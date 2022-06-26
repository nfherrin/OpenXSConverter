!OpenXSConverter is licensed under the MIT License.
!-------------------------------------------------------------------------------
! OpenXSConverter to convert different XS formats back and forth
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
PROGRAM openxsconverter
  USE globals
  USE infuncs
  USE outfuncs
  IMPLICIT NONE

  !read in command line arguments
  CALL readcl()
  WRITE(*,'(A)')'Reading in xs data'
  !read in xs file
  CALL readxs()
  IF(outformat .EQ. 'openmc')xsout=TRIM(xsout)//'.py'
  WRITE(*,'(2A)')'Outputting xs data to: ',TRIM(xsout)
  !output xs file
  CALL outputxs()

  WRITE(*,'(A)')'XS Converter completed! No detected errors in conversion.'
END PROGRAM openxsconverter
