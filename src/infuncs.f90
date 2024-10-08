!OpenXSConverter is licensed under the MIT License.
!input functions
MODULE infuncs
  USE HDF5
  USE globals
  USE string_module
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: readcl,readxs
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Gets command line arguments
  SUBROUTINE readcl()
    INTEGER::arg_count
    CHARACTER(64) :: tchar1
    arg_count=COMMAND_ARGUMENT_COUNT()

    IF(arg_count .GT. 3)STOP 'Only three argument variables are allowed for this program! &
        &"xsin" "outformat" "anis_ord_out"'

    !either use or prompt for input file name
    IF(arg_count .GE. 1)THEN
        CALL GET_COMMAND_ARGUMENT(1, xsin)
    ELSE
        WRITE(*,'(A)')'Input cross sections filename?'
        WRITE(*,'(A)',ADVANCE='NO')'> '
        READ(*,*)xsin
    END IF
    xsin=TRIM(ADJUSTL(xsin))

    !either use or prompt for output file name
    IF(arg_count .GE. 2)THEN
        CALL GET_COMMAND_ARGUMENT(2, outformat)
    ELSE
        WRITE(*,'(A)')'Output cross sections format? Available formats below:'
        WRITE(*,'(A)')'THOR'
        WRITE(*,'(A)')'OpenMC'
        WRITE(*,'(A)')'MCNP'
        WRITE(*,'(A)',ADVANCE='NO')'> '
        READ(*,*)outformat
    END IF
    outformat=TRIM(ADJUSTL(lowercase(outformat)))

    !either get the anisotropic scattering order or use what's in the file
    IF(arg_count .GE. 3)THEN
        CALL GET_COMMAND_ARGUMENT(3, tchar1)
        READ(tchar1,*)anis_out
    ELSE
        anis_out=-999
    END IF

    xsout=TRIM(xsin)//'_'//TRIM(outformat)//'.out'
  ENDSUBROUTINE readcl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Get xs format and call the xs reader based on format
  SUBROUTINE readxs()
    INTEGER :: ios,nwords,i
    CHARACTER(1000) :: tchar1,informat
    CHARACTER(100) :: words(200)

    !open xsin file
    xsin=TRIM(ADJUSTL(xsin))
    OPEN(UNIT=22,FILE=xsin,STATUS='OLD',ACTION='READ',IOSTAT=ios,IOMSG=tchar1)
    IF(ios .NE. 0)THEN
        WRITE(*,'(A)')tchar1
        STOP
    ENDIF

    !find the input format
    informat=''
    DO
      IF(xsin(LEN_TRIM(xsin)-2:LEN_TRIM(xsin)) .EQ. ".h5")THEN
        informat='openmc'
        EXIT
      ENDIF
      READ(22,*,IOSTAT=ios)tchar1
      SELECTCASE(tchar1)
        CASE('VERSION')
          BACKSPACE(22)
          READ(22,'(A20000)')tchar1
          CALL parse(tchar1,'=',words,nwords)
          tchar1=TRIM(ADJUSTL(words(2)))
          CALL parse(tchar1,"'",words,nwords)
          tchar1=TRIM(ADJUSTL(words(2)))
          CALL parse(tchar1," ",words,nwords)
          words(1)=TRIM(ADJUSTL(lowercase(words(1))))
          tchar1=TRIM(ADJUSTL(words(2)))
          IF(words(1) .EQ. 'serpent')THEN
            SELECTCASE(tchar1(1:1))
              CASE('1')
                informat='serp_gen_v1'
                EXIT
              CASE('2')
                informat='serp_gen_v2'
                EXIT
            ENDSELECT
          ENDIF
        CASE('THOR_XS_V1')
          informat='thor_v1'
          EXIT
        CASE DEFAULT
      ENDSELECT
      IF(ios .NE. 0)STOP 'No valid format identifier found'
    ENDDO

    !call the appropriate XS reader
    REWIND(22)
    SELECTCASE(informat)
      CASE('serp_gen_v1')
        CALL read_serp_v1()
      CASE('serp_gen_v2')
        CALL read_serp_v2()
      CASE('thor_v1')
        CALL read_thor_v1()
      CASE('openmc')
        CALL read_openmc()
      CASE  DEFAULT
        WRITE(*,'(3A)')'ERROR: ',TRIM(informat),' not a known input xs format.'
        STOP 'Fatal error'
    ENDSELECT
    IF(anis_out .LT. 0)THEN
      anis_out=levelanis
    ELSE
      IF(anis_out .GT. levelanis)THEN
        WRITE(*,'(A,I0,A,I0,A)')'WARNING: Requested Pn order ',anis_out,' greater than available Pn order ',levelanis,'!'
        WRITE(*,'(A,I0,A)')'Defaulting to maximum available Pn order ',levelanis,'.'
        anis_out=levelanis
      ENDIF
    ENDIF
    DO i=1,nummats
      IF(MAXVAL(chi(i,:)) .LE. 1.0D-12)chi(i,1)=1.0D0
      chi(i,:)=chi(i,:)/SUM(chi(i,:))
    ENDDO
    !close the input file
    CLOSE(22)
  ENDSUBROUTINE readxs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_serp_v1()
    STOP 'read_serp_v1 not complete'
  ENDSUBROUTINE read_serp_v1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_serp_v2()
    CHARACTER(100) :: tchar1
    INTEGER :: ios,m

    numgroups=0
    nummats=0
    levelanis=7
    ios=0
    !find the number of materials and number of energy groups
    DO WHILE(ios .EQ. 0)
      READ(22,*,IOSTAT=ios)tchar1
      !get the number of materials
      IF(tchar1 .EQ. "GC_UNIVERSE_NAME")nummats=nummats+1
      !get the number of energy groups
      IF(tchar1 .EQ. "MACRO_NG")THEN
        IF(numgroups .EQ. 0)THEN
          BACKSPACE(22)
          READ(22,*)tchar1,tchar1,tchar1,tchar1,numgroups
          ALLOCATE(eg_struc(numgroups+1))
          eg_struc=0.0
        ENDIF
      ENDIF
      !get the energy group structure
      IF(tchar1 .EQ. "MACRO_E")THEN
        IF(ALLOCATED(eg_struc))THEN
          BACKSPACE(22)
          READ(22,*)tchar1,tchar1,tchar1,tchar1,tchar1,tchar1,eg_struc(:)
        ENDIF
      ENDIF
    ENDDO
    REWIND(22)

    ALLOCATE(chi(nummats,numgroups),sigmaf(nummats,numgroups),nuf(nummats,numgroups))
    ALLOCATE(sigmat(nummats,numgroups),sigmas(nummats,levelanis+1,numgroups,numgroups))
    ALLOCATE(sigmaa(nummats,numgroups))

    m=0
    DO
      READ(22,*)tchar1
      !found the next material
      IF(tchar1 .EQ. 'GC_UNIVERSE_NAME')THEN
        m=m+1

        !absorption xs comes first
        CALL getserpv2xsdata('INF_ABS',sigmaa(m,:))
        !then fission xs
        CALL getserpv2xsdata('INF_FISS',sigmaf(m,:))
        !then nuf
        CALL getserpv2xsdata('INF_NUBAR',nuf(m,:))
        !then scattering production
        CALL getserpv2xsdata('INF_SCATTP0',sigmat(m,:))
        !then chi
        CALL getserpv2xsdata('INF_CHIT',chi(m,:))
        !then sigmass values for each scattering order
        CALL getserpv2scatdata('INF_SP0',sigmas(m,1,:,:))
        CALL getserpv2scatdata('INF_SP1',sigmas(m,2,:,:))
        CALL getserpv2scatdata('INF_SP2',sigmas(m,3,:,:))
        CALL getserpv2scatdata('INF_SP3',sigmas(m,4,:,:))
        CALL getserpv2scatdata('INF_SP4',sigmas(m,5,:,:))
        CALL getserpv2scatdata('INF_SP5',sigmas(m,6,:,:))
        CALL getserpv2scatdata('INF_SP6',sigmas(m,7,:,:))
        CALL getserpv2scatdata('INF_SP7',sigmas(m,8,:,:))
        !reset sigmat to not have an implicit negative abs XS
        sigmat(m,:)=sigmat(m,:)+sigmaa(m,:)
      ENDIF
      IF(m .GE. nummats)EXIT
    ENDDO
  ENDSUBROUTINE read_serp_v2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE getserpv2xsdata(xsID,xsvec)
    CHARACTER(*),INTENT(IN) :: xsID
    REAL(8) :: xsvec(numgroups)
    REAL(8) :: tvec(numgroups*2)
    CHARACTER(100) :: tchar1
    INTEGER :: g

    !get the xs data
    CALL find_line(xsID)
    READ(22,*)tchar1,tchar1,tchar1,tchar1,tchar1,tchar1,tvec(:)
    DO g=1,numgroups
      xsvec(g)=tvec(g*2-1)
    ENDDO
  ENDSUBROUTINE getserpv2xsdata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE getserpv2scatdata(xsID,xsvec)
    CHARACTER(*),INTENT(IN) :: xsID
    REAL(8) :: xsvec(numgroups,numgroups)
    REAL(8) :: tvec(2*numgroups**2)
    CHARACTER(100) :: tchar1
    INTEGER :: g,gp

    !get the xs data
    CALL find_line(xsID)
    READ(22,*)tchar1,tchar1,tchar1,tchar1,tchar1,tchar1,tvec
    DO g=1,numgroups
      DO gp=1,numgroups
        xsvec(gp,g)=tvec((g-1)*numgroups*2+gp*2-1)
      ENDDO
    ENDDO
  ENDSUBROUTINE getserpv2scatdata

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_thor_v1()
    INTEGER :: nwords,g,m,gp,l
    CHARACTER(100) :: words(2000)
    DO
      CALL get_thor_line(words,nwords)
      words(1)=TRIM(ADJUSTL(words(1)))
      IF(words(1) .EQ. 'THOR_XS_V1')THEN
        EXIT
      ENDIF
    ENDDO
    READ(words(2),*)nummats
    READ(words(3),*)numgroups
    READ(words(4),*)levelanis

    ALLOCATE(chi(nummats,numgroups),sigmaf(nummats,numgroups),nuf(nummats,numgroups))
    ALLOCATE(sigmat(nummats,numgroups),sigmas(nummats,levelanis+1,numgroups,numgroups))
    ALLOCATE(eg_struc(numgroups),sigmaa(nummats,numgroups))
    eg_struc=0.0

    CALL get_thor_line(words,nwords)
    words(1)=TRIM(ADJUSTL(lowercase(words(1))))
    IF(words(1) .NE. 'id')THEN
      DO g=1,numgroups
        READ(words(g),*)eg_struc(g)
      ENDDO
    ELSE
      !generate logarithmically even groups
      DO m=1,numgroups
        eg_struc(m)=10**(8.0+(m-1.0)*(-14.0D0)/(numgroups-1))
      ENDDO
    ENDIF
    REWIND(22)
    !get all materials data
    DO m=1,nummats
      DO
        !get next line
        CALL get_thor_line(words,nwords)
        words(1)=TRIM(ADJUSTL(lowercase(words(1))))
        IF(words(1) .EQ. 'id')THEN
          !read in the fission spectrum
          CALL get_thor_line(words,nwords)
          DO g=1,numgroups
            READ(words(g),*)chi(m,g)
          ENDDO
          !read in SigmaF
          CALL get_thor_line(words,nwords)
          DO g=1,numgroups
            READ(words(g),*)sigmaf(m,g)
          ENDDO
          !read in nu
          CALL get_thor_line(words,nwords)
          DO g=1,numgroups
            READ(words(g),*)nuf(m,g)
          ENDDO
          !read in total/transport xs
          CALL get_thor_line(words,nwords)
          DO g=1,numgroups
            READ(words(g),*)sigmat(m,g)
          ENDDO
          !read in scattering xs format gp->g
          DO l=1,levelanis+1
            DO g=1,numgroups
              CALL get_thor_line(words,nwords)
              DO gp=1,numgroups
                READ(words(gp),*)sigmas(m,l,g,gp)
              ENDDO
            ENDDO
          ENDDO
          EXIT
        ENDIF
      ENDDO
    ENDDO

    CALL comp_sig_a()
  ENDSUBROUTINE read_thor_v1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_openmc()
    INTEGER :: err1,m
    INTEGER(HID_T) :: file_id,dataset_id,dataspace_id,temp_hid1(3),temp_hid2(3)
    CHARACTER(64) :: mc_group,xs_dtype

    CALL h5open_f(err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error initializing HDF5 routines"
    ENDIF

    CALL h5fopen_f(xsin, H5F_ACC_RDONLY_F, file_id, err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening HDF5 file"
    ENDIF

    CALL h5gget_obj_info_idx_f(file_id,'/',0,xs_dtype,m,err1)
    xs_dtype=TRIM(ADJUSTL(xs_dtype))
    IF (err1 .LT. 0) THEN
      STOP " *** Error getting XS datatype"
    ENDIF

    !determine number of materials/cells that are being indexed
    CALL h5gn_members_f(file_id,'/'//TRIM(xs_dtype)//'/',nummats,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error reading in materials/cells numbers"
    ELSEIF(nummats .LE. 0)THEN
      STOP ' *** Error less that 1 materials/cells found'
    ENDIF

    !determine number of energy groups
    CALL h5dopen_f(file_id,'/'//TRIM(xs_dtype)//'/1/total/average',dataset_id,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening HDF5 dataset"
    ENDIF
    CALL h5dget_space_f(dataset_id,dataspace_id,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening HDF5 dataspace"
    ENDIF
    CALL h5sget_simple_extent_ndims_f(dataspace_id,m,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening getting level of anisotropy 1"
    ENDIF
    CALL h5sget_simple_extent_dims_f(dataspace_id,temp_hid1,temp_hid2,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error getting number of energy groups"
    ENDIF
    numgroups=INT(temp_hid1(1),4)

    !determine level of anisotropy
    CALL h5dopen_f(file_id,'/'//TRIM(xs_dtype)//'/1/scatter matrix/average',dataset_id,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening HDF5 dataset"
    ENDIF
    CALL h5dget_space_f(dataset_id,dataspace_id,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening getting HDF5 dataspace"
    ENDIF
    CALL h5sget_simple_extent_ndims_f(dataspace_id,m,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening getting level of anisotropy"
    ENDIF
    IF(m .EQ. 2)THEN
      levelanis=0
    ELSEIF(m .EQ. 3)THEN
      CALL h5sget_simple_extent_dims_f(dataspace_id,temp_hid1,temp_hid2,err1)
      IF (err1 .LT. 0) THEN
        STOP " *** Error getting number of energy groups"
      ENDIF
      levelanis=INT(temp_hid1(1)-1,4)
    ELSE
      STOP 'bad number of anisotropy determinants, only use legendre orders!'
    ENDIF

    ALLOCATE(chi(nummats,numgroups),sigmaf(nummats,numgroups),nuf(nummats,numgroups))
    ALLOCATE(sigmat(nummats,numgroups),sigmas(nummats,levelanis+1,numgroups,numgroups))
    ALLOCATE(eg_struc(numgroups),sigmaa(nummats,numgroups))
    eg_struc=0.0

    !generate logarithmically even groups (we store in terms of MeV)
    WRITE(*,'(A)')'WARNING: Energy data is not read in conversion from OpenMC'
    WRITE(*,'(A)')'Rather it is logarithmically generated from 20 MeV to 0.02 eV.'
    WRITE(*,'(A)')'User should adjust output accordingly if energy group structure is important for the application.'
    DO m=1,numgroups
      eg_struc(m)=2.0D0*10**(1.0D0+(m-1.0D0)*(-9.0D0)/(numgroups-1.0D0))
    ENDDO

    DO m=1,nummats
      WRITE(mc_group,'(A,I0,A)')'/'//TRIM(xs_dtype)//'/',m,'/chi/average'
      CALL readin_omc_xs(file_id,TRIM(mc_group),chi(m,:))
      WRITE(mc_group,'(A,I0,A)')'/'//TRIM(xs_dtype)//'/',m,'/fission/average'
      CALL readin_omc_xs(file_id,TRIM(mc_group),sigmaf(m,:))
      WRITE(mc_group,'(A,I0,A)')'/'//TRIM(xs_dtype)//'/',m,'/nu-fission/average'
      CALL readin_omc_xs(file_id,TRIM(mc_group),nuf(m,:))
      IF(MINVAL(sigmaf(m,:)) .GT. 0.0D0)nuf(m,:)=nuf(m,:)/sigmaf(m,:)
      WRITE(mc_group,'(A,I0,A)')'/'//TRIM(xs_dtype)//'/',m,'/total/average'
      CALL readin_omc_xs(file_id,TRIM(mc_group),sigmat(m,:))
      WRITE(mc_group,'(A,I0,A)')'/'//TRIM(xs_dtype)//'/',m,'/scatter matrix/average'
      CALL readin_omc_smat(file_id,TRIM(mc_group),sigmas(m,:,:,:))
    ENDDO

    CALL comp_sig_a()
  ENDSUBROUTINE read_openmc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE readin_omc_xs(file_id,group_char,xs_arr)
    INTEGER(HID_T),INTENT(IN) :: file_id
    CHARACTER(*),INTENT(IN) :: group_char
    REAL(8),INTENT(OUT) :: xs_arr(:)
    INTEGER(HID_T) :: dataset_id,datatype_id,dataspace_id
    INTEGER(HSIZE_T) :: dims(1)
    INTEGER :: err1

    CALL h5dopen_f(file_id,group_char,dataset_id,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening HDF5 dataset"
    ENDIF
    CALL h5dget_space_f(dataset_id,dataspace_id,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening HDF5 dataspace"
    ENDIF
    CALL h5dget_type_f(dataset_id,datatype_id,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening HDF5 dataspace"
    ENDIF

    xs_arr=0.0
    dims=numgroups
    CALL h5dread_f(dataset_id,datatype_id,xs_arr,dims,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening HDF5 data"
    ENDIF
  ENDSUBROUTINE readin_omc_xs

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE readin_omc_smat(file_id,group_char,scat_arr)
    INTEGER(HID_T),INTENT(IN) :: file_id
    CHARACTER(*),INTENT(IN) :: group_char
    REAL(8),INTENT(OUT) :: scat_arr(:,:,:)
    INTEGER(HID_T) :: dataset_id,datatype_id,dataspace_id
    INTEGER(HSIZE_T) :: dims(3)
    INTEGER :: err1

    CALL h5dopen_f(file_id,group_char,dataset_id,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening HDF5 dataset"
    ENDIF
    CALL h5dget_space_f(dataset_id,dataspace_id,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening HDF5 dataspace"
    ENDIF
    CALL h5dget_type_f(dataset_id,datatype_id,err1)
    IF (err1 .LT. 0) THEN
      STOP " *** Error opening HDF5 dataspace"
    ENDIF

    scat_arr=0.0
    dims=numgroups
    IF(levelanis .EQ. 0)THEN
      CALL h5dread_f(dataset_id,datatype_id,scat_arr(1,:,:),dims(1:2),err1)
      IF (err1 .LT. 0) THEN
        STOP " *** Error opening HDF5 data"
      ENDIF
    ELSE
      dims(1)=levelanis+1
      CALL h5dread_f(dataset_id,datatype_id,scat_arr,dims,err1)
      IF (err1 .LT. 0) THEN
        STOP " *** Error opening HDF5 data"
      ENDIF
    ENDIF
  ENDSUBROUTINE readin_omc_smat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_thor_line(words,nwords)
    INTEGER,INTENT(OUT) :: nwords
    CHARACTER(100),INTENT(OUT) :: words(2000)
    CHARACTER(20000) :: twords(3)
    CHARACTER(20000) :: tchar1
    INTEGER :: ios
    DO
      READ(22,'(A10000)',IOSTAT=ios)tchar1
      IF( ios .NE. 0)STOP 'end of xs file was reached before all data/materials were found'
      tchar1=TRIM(ADJUSTL(tchar1))
      !finding uncommented line that isn't empty
      IF(tchar1(1:1) .NE. '!' .AND. tchar1 .NE. '')THEN
        !ignore commented portions of line
        CALL parse(tchar1,'!',twords,nwords)
        tchar1=TRIM(ADJUSTL(twords(1)))
        CALL parse(tchar1,' ',words,nwords)
        EXIT
      ENDIF
    ENDDO
  ENDSUBROUTINE get_thor_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE find_line(linechar)
    CHARACTER(*),INTENT(IN) :: linechar
    CHARACTER(80) :: tchar1

    DO
      READ(22,*)tchar1
      IF(TRIM(ADJUSTL(tchar1)) .EQ. TRIM(ADJUSTL(linechar)))EXIT
    ENDDO
    BACKSPACE(22)
  ENDSUBROUTINE find_line

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE comp_sig_a()
    INTEGER :: m,g,gp

    sigmaa=sigmat
    DO m=1,nummats
      DO g=1,numgroups
        DO gp=1,numgroups
          sigmaa(m,gp)=sigmaa(m,gp)-sigmas(m,1,g,gp)
        ENDDO
      ENDDO
    ENDDO
  ENDSUBROUTINE comp_sig_a
END MODULE infuncs
