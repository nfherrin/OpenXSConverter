!OpenXSConverter is licensed under the MIT License.
!output functions
MODULE outfuncs
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: outputxs
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !based on xs output format, call the xs out-putter
  SUBROUTINE outputxs()
    INTEGER :: i

    DO i=1,numgroups
      IF(eg_struc(i) .GE. 20.0D0)eg_struc(i)=20.0D0
    ENDDO
    SELECTCASE(outformat)
      CASE('thor')
        CALL out_thor()
      CASE('mcnp')
        CALL out_mcnp()
      CASE('openmc')
        CALL out_openmc()
      CASE('moose')
        CALL out_moose()
      CASE DEFAULT
        STOP 'bad output format'
    ENDSELECT
  ENDSUBROUTINE outputxs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE out_thor()
    INTEGER :: ios,m,gp,l
    CHARACTER(10000) :: tchar1

    !open xsout file
    OPEN(UNIT=32,FILE=xsout,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ios,IOMSG=tchar1)
    IF(ios .NE. 0)THEN
        WRITE(*,'(A)')tchar1
        STOP
    ENDIF

    !output the xs characteristics
    WRITE(32,'(A,I0,A,I0,A,I0)')'THOR_XS_V1 ',nummats,' ',numgroups,' ',anis_out
    !output the energy group structure
    WRITE(tchar1,'(10000ES20.12)')eg_struc(1:numgroups)
    WRITE(32,'(A,ES20.12)')TRIM(ADJUSTL(tchar1)),0.0
    !output the xs data
    DO m=1,nummats
      WRITE(32,'(A,I0)')'id ',m
      WRITE(tchar1,'(10000ES16.8)')chi(m,:)
      WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
      WRITE(tchar1,'(10000ES16.8)')sigmaf(m,:)
      WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
      WRITE(tchar1,'(10000ES16.8)')nuf(m,:)
      WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
      WRITE(tchar1,'(10000ES16.8)')sigmat(m,:)
      WRITE(32,'(A)')TRIM(ADJUSTL(tchar1))
      DO l=1,anis_out+1
        DO gp=1,numgroups
          WRITE(32,'(10000ES16.8)')sigmas(m,l,gp,:)
        ENDDO
      ENDDO
    ENDDO
    !close the output file
    CLOSE(32)
  ENDSUBROUTINE out_thor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!mcnp output, see https://mcnp.lanl.gov/pdf_files/la-12704.pdf for the format of the multigroup transport table
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE out_mcnp()
    CHARACTER(64)::tempcharacter
    REAL(8),ALLOCATABLE::xsarray(:),equi_bins(:,:,:,:)
    INTEGER::ios,i,j,g,l,arrloc
    !NXS values
    INTEGER :: LDB,NLEG
    !JXS values
    INTEGER :: LERG,LTOT,LFISS,LNU,LCHI,LABS,LP0L,LXPNL,LPNL

    !20 equi-probable cosine bins, change this if you want more/less
    NLEG=21
    ALLOCATE(equi_bins(nummats,numgroups,numgroups,NLEG))
    equi_bins=0.0D0
    CALL compute_equi_cos_bins(equi_bins,NLEG)
    LFISS=0
    LNU=0
    LCHI=0

    DO i=1,nummats
      WRITE(xsout,'(A,I0)')'xs_'//TRIM(xsin)//'_'//TRIM(outformat)//'.mat',i
      !open xsout file
      OPEN(UNIT=32,FILE=xsout,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ios,IOMSG=tempcharacter)
      IF(ios .NE. 0)THEN
        WRITE(*,*)tempcharacter
        STOP
      END IF

      LDB=(2+NLEG)*numgroups**2+numgroups*7+3

      !build xs array (very complex, see MCNP5 manual Appendix F for details)
      ALLOCATE(xsarray(LDB))
      xsarray=0
      arrloc=1

      !first block is the energy bounds
      LERG=arrloc
      DO j=1,numgroups
        xsarray(arrloc)=(eg_struc(j)+eg_struc(j+1))/2.0D0
        arrloc=arrloc+1
      ENDDO

      !second block is the energy group widths
      DO j=1,numgroups
        xsarray(arrloc)=eg_struc(j)-eg_struc(j+1)
        arrloc=arrloc+1
      ENDDO

      !third block is the total cross section
      LTOT=arrloc
      DO j=1,numgroups
          xsarray(arrloc)=sigmat(i,j)
          arrloc=arrloc+1
      ENDDO

      !fourth block is fission cross section
      LFISS=arrloc
      DO j=1,numgroups
        xsarray(arrloc)=sigmaf(i,j)
        arrloc=arrloc+1
      ENDDO

      !fifth block is nu
      LNU=arrloc
      DO j=1,numgroups
        xsarray(arrloc)=nuf(i,j)
        arrloc=arrloc+1
      ENDDO

      !sixth block is fission spectrum
      LCHI=arrloc
      DO j=1,numgroups
        xsarray(arrloc)=chi(i,j)
        arrloc=arrloc+1
      ENDDO

      !seventh block is absorption XS
      LABS=arrloc
      DO j=1,numgroups
        xsarray(arrloc)=sigmaa(i,j)
        arrloc=arrloc+1
      ENDDO

      !eighth block is P0 scattering
      LP0L=arrloc
      xsarray(arrloc)=arrloc+1
      arrloc=arrloc+1
      DO j=1,numgroups
        DO g=1,numgroups
          xsarray(arrloc)=sigmas(i,1,g,j)
          arrloc=arrloc+1
        ENDDO
      ENDDO

      !ninth block is the XPN block
      LXPNL=arrloc
      xsarray(arrloc)=arrloc+1
      arrloc=arrloc+1
      ios=1
      DO j=1,numgroups
        DO g=1,numgroups
          xsarray(arrloc)=ios
          ios=ios+NLEG
          arrloc=arrloc+1
        ENDDO
      ENDDO

      !tenth block is the PN scattering
      LPNL=arrloc
      xsarray(arrloc)=arrloc+1
      arrloc=arrloc+1
      DO j=1,numgroups
        DO g=1,numgroups
          DO l=1,NLEG
            xsarray(arrloc)=equi_bins(i,g,j,l)
            arrloc=arrloc+1
          ENDDO
        ENDDO
      ENDDO

      !print out informative data
      WRITE(32,'(A,I0,A)')'  ',1110+i,'.00m  1.0 '
      WRITE(32,'(A,I0,A,I0,3A,I0,A,ES11.4)')'  in MCNP: xs',i,' 111',i,'.00m 1.0 ',TRIM(xsout)&
          &,' 0 1 1 ',LDB,' 0 0',0.0
      WRITE(32,*)
      WRITE(32,*)
      WRITE(32,*)
      WRITE(32,*)
      !output the NXS array
      WRITE(32,'(8I9)')LDB,1110+i,NLEG,0,numgroups,numgroups-1,numgroups-1,0
      WRITE(32,'(8I9)')0,1,0,1,0,0,0,0
      !output the JXS array
      WRITE(32,'(8I9)')LERG,LTOT,LFISS,LNU,LCHI,LABS,0,0
      WRITE(32,'(8I9)')0,0,0,0,LP0L,0,0,LXPNL
      WRITE(32,'(8I9)')LPNL,0,0,0,0,0,0,0
      WRITE(32,'(8I9)')0,0,0,0,0,0,0,0
      WRITE(32,'(4ES20.12)',ADVANCE='NO')xsarray(:)
      WRITE(32,*)
      DEALLOCATE(xsarray)

      CLOSE(32)
    ENDDO

  ENDSUBROUTINE out_mcnp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!compute equi-probable cosine bin boundaries based on the anisotropic approximation given
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE compute_equi_cos_bins(equi_bins,NLEG)
    REAL(8), INTENT(OUT) :: equi_bins(:,:,:,:)
    INTEGER, INTENT(IN) :: NLEG
    REAL(8) :: h,x,x_old,int_1,int_2
    INTEGER :: i,g,j,l,count

    h=1.0D-5 !this is what we will step integrate with using simpson's rule until we reach the next segment
    equi_bins=0.0D0
    DO i=1,nummats
      DO g=1,numgroups
        DO j=1,numgroups
          IF(sigmas(i,1,g,j) .GT. 0.0D0 .AND. anis_out .NE. 0)THEN
            !generates the probability cosine bins based on the given cross sections for the left half
            equi_bins(i,g,j,1)=-1.0D0
            x_old=-1.0D0
            x=x_old
            DO l=2,(NLEG-1)/2
              count=0
              int_1=0.0D0
              DO
                int_1=int_1+h/6.0D0*(xs_normalized(i,g,j,x)+xs_normalized(i,g,j,x+h)+4.0*xs_normalized(i,g,j,x+h/2.0D0))
                count=count+1
                x=x_old+count*h
                IF(1.0D0/(NLEG-1.0D0)-int_1 .LT. 0.0D0)THEN
                  int_2=int_1-h/6.0D0*(xs_normalized(i,g,j,x-h)+xs_normalized(i,g,j,x)+4.0*xs_normalized(i,g,j,x-h/2.0D0))
                  !interpolate to get a better result
                  x=x-(1.0D0/(NLEG-1.0D0)-int_1)*h/(int_2-int_1)
                  EXIT
                ENDIF
              ENDDO
              x_old=x
              equi_bins(i,g,j,l)=x
            ENDDO
            !generates the probability cosine bins based on the given cross sections for the right half (splitting like this helps avoid accumulation of roundoff error)
            equi_bins(i,g,j,NLEG)=1.0D0
            x_old=1.0D0
            x=x_old
            DO l=NLEG-1,(NLEG-1)/2+1,-1
              count=0
              int_1=0.0D0
              DO
                int_1=int_1+h/6.0D0*(xs_normalized(i,g,j,x)+xs_normalized(i,g,j,x-h)+4.0*xs_normalized(i,g,j,x-h/2.0D0))
                count=count+1
                x=x_old-count*h
                IF(1.0D0/(NLEG-1.0D0)-int_1 .LT. 0.0D0)THEN
                  int_2=int_1-h/6.0D0*(xs_normalized(i,g,j,x+h)+xs_normalized(i,g,j,x)+4.0*xs_normalized(i,g,j,x+h/2.0D0))
                  x=x+(1.0D0/(NLEG-1.0D0)-int_1)*h/(int_2-int_1)
                  EXIT
                ENDIF
              ENDDO
              x_old=x
              equi_bins(i,g,j,l)=x
            ENDDO
          ELSE
            !if the scattering is 0, set it all to equal
            DO l=1,NLEG
              equi_bins(i,g,j,l)=-1.0D0+(l-1.0D0)*2.0D0/(NLEG-1.0D0)
            ENDDO
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDSUBROUTINE compute_equi_cos_bins

  !computes a Legendre polynomial using the recurrence relation
  REAL(8) FUNCTION p_l(n,x)
    REAL(8), INTENT(IN) :: x
    INTEGER, INTENT(IN) :: n
    REAL(8) :: p_l_n(n+1)
    INTEGER :: i

    p_l_n(1)=1.0D0
    IF(n .GT. 0)THEN
      p_l_n(2)=x
      DO i=1,n-1
        p_l_n(i+2)=((2.0D0*i+1.0D0)*x*p_l_n(i+1) - i*p_l_n(i))/(i+1.0D0)
      ENDDO
    ENDIF
    p_l=p_l_n(n+1)
  ENDFUNCTION p_l

  !computes the normalized value of the cross section at point x
  REAL(8) FUNCTION xs_normalized(i,g,j,x)
    REAL(8),INTENT(IN) :: x
    INTEGER, INTENT(IN) :: i,g,j
    INTEGER :: l

    xs_normalized=0.0D0
    DO l=1,anis_out+1
      xs_normalized=xs_normalized+sigmas(i,l,g,j)*p_l(l-1,x)/(sigmas(i,1,g,j)*2.0D0)
    ENDDO
  ENDFUNCTION xs_normalized

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE out_openmc()
    INTEGER :: ios,g,m,gp,l
    CHARACTER(64) :: tchar1

    !open xsout file
    OPEN(UNIT=32,FILE=xsout,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ios,IOMSG=tchar1)
    IF(ios .NE. 0)THEN
        WRITE(*,'(A)')tchar1
        STOP
    ENDIF
    WRITE(32,'(A)')'import numpy as np'
    WRITE(32,'(A)')'import openmc'
    WRITE(32,'(A)')'import os'
    WRITE(32,'(A)')
    WRITE(32,'(A)')'#INSTRUCTIONS: place these cross sections at the top of your OpenMC python script'
    WRITE(32,'(A)')'#OR build your script from this baseline'
    WRITE(32,'(A)')'####################################################################################################'
    WRITE(32,'(A)')'###################################-Beginning of Cross Sections-####################################'
    WRITE(32,'(A)')'####################################################################################################'
    WRITE(32,'(A)')'#Group structure:'
    WRITE(32,'(A,ES15.8)',ADVANCE='NO')'groups = openmc.mgxs.EnergyGroups([',0.0D0
    DO g=0,numgroups-1
      IF(MOD(g+1,6) .EQ. 0)THEN
        WRITE(32,'(A)')''
        WRITE(32,'(A)',ADVANCE='NO')'                          '
      ENDIF
      !openmc input is in eV but we store in MeV, so convert
      WRITE(32,'(A,ES15.8)',ADVANCE='NO')',',eg_struc(numgroups-g)*10**6
    ENDDO
    WRITE(32,'(A)')'])'
    DO m=1,nummats
      WRITE(32,'(A,I0,A)')'#Data for Material ',m,':'
      WRITE(32,'(A,I0,A,I0,A)')"mat",m,"_xsdat = openmc.XSdata('mat_",m,"', groups)"
      WRITE(32,'(A,I0,A,I0)')"mat",m,"_xsdat.order = ",anis_out
      !total xs
      CALL print_xs_openmc(m,'total',sigmat(m,:))
      CALL print_xs_openmc(m,'absorption',sigmaa(m,:))
      CALL print_xs_openmc(m,'fission',sigmaf(m,:))
      CALL print_xs_openmc(m,'nu_fission',nuf(m,:)*sigmaf(m,:))
      CALL print_xs_openmc(m,'chi',chi(m,:))
      !print the scattering matrix, a bit more involved...
      WRITE(32,'(A)')'scatter_matrix = np.array(\'
      WRITE(32,'(A)',ADVANCE='NO')'    ['
      DO l=1,anis_out+1
        IF(l .NE. 1)WRITE(32,'(A)',ADVANCE='NO')'     '
        WRITE(32,'(A)',ADVANCE='NO')'['
        DO gp=1,numgroups
          IF(gp .NE. 1)WRITE(32,'(A)',ADVANCE='NO')'      '
          WRITE(32,'(A)',ADVANCE='NO')'['
          DO g=1,numgroups
            IF(g .NE. 1)WRITE(32,'(A)',ADVANCE='NO')','
            WRITE(32,'(ES15.8)',ADVANCE='NO')sigmas(m,l,gp,g)
          ENDDO
          WRITE(32,'(A)',ADVANCE='NO')']'
          IF(gp .NE. numgroups)WRITE(32,'(A)')','
        ENDDO
        WRITE(32,'(A)',ADVANCE='NO')']'
        IF(l .NE. anis_out+1)WRITE(32,'(A)')','
      ENDDO
      WRITE(32,'(A)')'])'
      WRITE(32,'(A)')'scatter_matrix = np.transpose(scatter_matrix)'
      WRITE(32,'(A,I0,A)')"mat",m,"_xsdat.set_scatter_matrix(scatter_matrix, temperature=294.)"
    ENDDO
    WRITE(32,'(A)')'#Create the cross sections hdf5 file:'
    WRITE(32,'(A)')'mg_cross_sections_file = openmc.MGXSLibrary(groups)'
    DO m=1,nummats
      WRITE(32,'(A,I0,A)')'mg_cross_sections_file.add_xsdata(mat',m,'_xsdat)'
    ENDDO
    WRITE(32,'(A)')"mg_cross_sections_file.export_to_hdf5('cross_sections.h5')"
    WRITE(32,'(A)')'#Assign each cross section to a separate material:'
    WRITE(32,'(A)')'materials = {}'
    WRITE(32,'(A)',ADVANCE='NO')"for xs in ["
    DO m=1,nummats
      IF(m .NE. 1)WRITE(32,'(A)',ADVANCE='NO')','
      WRITE(32,'(A,I0,A)',ADVANCE='NO')"'mat_",m,"'"
    ENDDO
    WRITE(32,'(A)')']:'
    WRITE(32,'(A)')"    materials[xs] = openmc.Material(name=xs)"
    WRITE(32,'(A)')"    materials[xs].set_density('macro', 1.)"
    WRITE(32,'(A)')"    materials[xs].add_macroscopic(xs)"
    WRITE(32,'(A)')'#Create the materials file for this specification:'
    WRITE(32,'(A)')"materials_file = openmc.Materials(materials.values())"
    WRITE(32,'(A)')"materials_file.cross_sections = 'cross_sections.h5'"
    WRITE(32,'(A)')"materials_file.export_to_xml()"
    WRITE(32,'(A)')'####################################################################################################'
    WRITE(32,'(A)')'######################################-End of Cross Sections-#######################################'
    WRITE(32,'(A)')'####################################################################################################'
    WRITE(32,'(A)')"#The following creates the settings and sets the energy mode to multi-group:"
    WRITE(32,'(A)')"settings_file = openmc.Settings()"
    WRITE(32,'(A)')"settings_file.energy_mode = 'multi-group'"
    WRITE(32,'(A)')"#The user should complete the OpenMC input below by specifying geometry, settings, etc."
  ENDSUBROUTINE out_openmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE print_xs_openmc(matnum,xstype,xsarr)
    INTEGER,INTENT(IN) :: matnum
    CHARACTER(*),INTENT(IN) :: xstype
    REAL(8),INTENT(IN) :: xsarr(*)
    INTEGER :: g

    WRITE(32,'(A,I0,3A,ES15.8)',ADVANCE='NO')'mat',matnum,'_xsdat.set_',TRIM(ADJUSTL(xstype)),'([',xsarr(1)
    DO g=2,numgroups
      IF(MOD(g,6) .EQ. 0)THEN
        WRITE(32,'(A)')''
        WRITE(32,'(A)',ADVANCE='NO')'                          '
      ENDIF
      WRITE(32,'(A,ES15.8)',ADVANCE='NO')',',xsarr(g)
    ENDDO
    WRITE(32,'(A)')'], temperature=294.)'
  ENDSUBROUTINE print_xs_openmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE out_moose()
    INTEGER :: ios,m,g
    CHARACTER(64) :: tchar1

    WRITE(*,'(A)')'NOTE: Given diffusion cross sections for MOOSE output are 1/(3*Sigma_t)'
    IF(anis_out .GT. 0)THEN
      WRITE(*,'(A,I0,A)')'WARNING! ',anis_out,' specified anisotropic XS, but for now MOOSE output &
        &only supports isotropic scattering.'
    ENDIF

    !open xsout file
    OPEN(UNIT=32,FILE=xsout,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ios,IOMSG=tchar1)
    IF(ios .NE. 0)THEN
        WRITE(*,'(A)')tchar1
        STOP
    ENDIF

    WRITE(32,'(A)')'<Materials>'
    WRITE(32,'(A,I0,A)')'   <Macros NG="',numgroups,'">'
    DO m=1,nummats
      IF(MAXVAL(ABS(sigmaf(m,:))) .GT. 0)THEN
        WRITE(32,'(A,I0,A)')'     <material ID="',m,'" fissile="true">'
      ELSE
        WRITE(32,'(A,I0,A)')'     <material ID="',m,'" fissile="false">'
      ENDIF
      WRITE(32,'(A,I0,A)')'        <name>mat_',m,'</name>'
      WRITE(32,'(A)')'        <TotalXS>'
      DO g=1,numgroups
        WRITE(32,'(A,ES16.8)')'        ',sigmat(m,g)
      ENDDO
      WRITE(32,'(A)')'        </TotalXS>'
      IF(MAXVAL(ABS(sigmaf(m,:))) .GT. 0)THEN
        WRITE(32,'(A)')'        <NuFissionXS>'
        DO g=1,numgroups
          WRITE(32,'(A,ES16.8)')'        ',nuf(m,g)*sigmaf(m,g)
        ENDDO
        WRITE(32,'(A)')'        </NuFissionXS>'
        WRITE(32,'(A)')'        <ChiXS>'
        DO g=1,numgroups
          WRITE(32,'(A,ES16.8)')'        ',chi(m,g)
        ENDDO
        WRITE(32,'(A)')'        </ChiXS>'
        WRITE(32,'(A)')'        <FissionXS>'
        DO g=1,numgroups
          WRITE(32,'(A,ES16.8)')'        ',sigmaf(m,g)
        ENDDO
        WRITE(32,'(A)')'        </FissionXS>'
      ENDIF
      WRITE(32,'(A)')'        <ScatteringXS>'
      DO g=1,numgroups
        WRITE(32,'(A,10000ES16.8)')'        ',sigmas(m,1,g,:)
      ENDDO
      WRITE(32,'(A)')'        </ScatteringXS>'
      WRITE(32,'(A)')'        <DiffusionCoefficient>'
      DO g=1,numgroups
        WRITE(32,'(A,ES16.8)')'        ',1.0D0/(3.0D0*sigmat(m,g))
      ENDDO
      WRITE(32,'(A)')'        </DiffusionCoefficient>'
      WRITE(32,'(A)')'     </material>'
    ENDDO
    WRITE(32,'(A)')'    </Macros>'
    WRITE(32,'(A)')'</Materials>'

    !close the output file
    CLOSE(32)
  ENDSUBROUTINE out_moose
END MODULE outfuncs
