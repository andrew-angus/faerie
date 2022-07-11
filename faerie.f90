MODULE faerie

USE ISO_C_BINDING
USE ISO_FORTRAN_ENV
USE NETCDF

IMPLICIT NONE

! Parameters
INTEGER, PARAMETER :: cd = C_DOUBLE
INTEGER, PARAMETER :: ci = C_INT
INTEGER, PARAMETER :: cb = C_BOOL
INTEGER, PARAMETER :: dp = REAL64
INTEGER, PARAMETER :: fi = INT32
REAL(dp), PARAMETER :: pi = ACOS(-1.0_dp)
REAL(dp), PARAMETER :: small = 1e-30_dp
REAL(dp), PARAMETER :: releps = 1.0e-6_dp
CHARACTER(LEN=1), PARAMETER :: creturn = ACHAR(13)

! Kernel function interface
INTERFACE
  FUNCTION kernel()
    USE ISO_FORTRAN_ENV
    REAL(REAL64) :: kernel
  END FUNCTION
END INTERFACE

! GP object
! single y output, zero prior mean, for now
TYPE gp_obj
  REAL(dp) :: noise, var
  REAL(dp), DIMENSION(:), ALLOCATABLE :: l
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: x
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: y
  INTEGER(fi) :: nx
  PROCEDURE(kernel), POINTER, NOPASS :: kern => NULL()
  !PROCEDURE(xconvert), POINTER, NOPASS :: xcon => NULL()
  !PROCEDURE(yconvert), POINTER, NOPASS :: ycon => NULL()
  !PROCEDURE(xrevert), POINTER, NOPASS :: xrev => NULL()
  !PROCEDURE(yrevert), POINTER, NOPASS :: yrev => NULL()
END TYPE
TYPE(gp_obj) :: gp

! Shared data
INTEGER(fi) :: i,j
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: dist_bnds

CONTAINS

! Read in dataset, distributions, and gp params
SUBROUTINE read_data()

  ! Known labels
  CHARACTER (LEN = *), PARAMETER :: fname = "gp.nc"
  CHARACTER (LEN = *), DIMENSION(4), PARAMETER :: &
    dimlabs = (/"inputs ","outputs","samples","bounds "/)
  CHARACTER (LEN = *), DIMENSION(4), PARAMETER :: &
    varlabs = (/"lengthscales  ","input_samples ","output_samples", &
                "dist_bounds   "/)

  ! IDs and lengths
  INTEGER :: fid,attlen
  INTEGER, DIMENSION(4) :: dimids, dimlens, varids

  ! Local output
  CHARACTER(:),ALLOCATABLE :: kern

  ! Open file
  CALL check(NF90_OPEN(fname, NF90_NOWRITE, fid))

  ! Get global attributes
  CALL check(NF90_INQUIRE_ATTRIBUTE(fid,NF90_GLOBAL,"kernel",len=attlen))
  ALLOCATE(CHARACTER(attlen)::kern)
  CALL check(NF90_GET_ATT(fid,NF90_GLOBAL,"kernel",kern))
  CALL check(NF90_GET_ATT(fid,NF90_GLOBAL,"var",gp%var))
  CALL check(NF90_GET_ATT(fid,NF90_GLOBAL,"noise",gp%noise))

  ! Get dimension/variable IDs and lengths
  DO i = 1, 4
    CALL check(NF90_INQ_DIMID(fid,TRIM(dimlabs(i)),dimids(i)))
    CALL check(NF90_INQUIRE_DIMENSION(fid,dimids(i),len=dimlens(i)))
    CALL check(NF90_INQ_VARID(fid,TRIM(varlabs(i)),varids(i)))
  END DO

  ! Read variables
  ALLOCATE(gp%l(dimlens(1)),gp%x(dimlens(1),dimlens(3)),gp%y(dimlens(2),dimlens(3)))
  ALLOCATE(dist_bnds(dimlens(4),dimlens(1)))
  CALL check(NF90_GET_VAR(fid, varids(1), gp%l))
  CALL check(NF90_GET_VAR(fid, varids(2), gp%x))
  CALL check(NF90_GET_VAR(fid, varids(3), gp%y))
  CALL check(NF90_GET_VAR(fid, varids(4), dist_bnds))

  ! Close the file, freeing all resources.
  CALL check(NF90_CLOSE(fid))

  ! Assign GP function pointer
  IF (kern == "rbf") THEN
    gp%kern => rbf
  ELSE IF (kern == "Mat52") THEN
    gp%kern => matern52
  ELSE IF (kern == "Mat32") THEN
    gp%kern => matern32
  ELSE IF (kern == "Exponential") THEN
    gp%kern => exponential
  ELSE
    STOP "Read kernel name invalid: must be one of rbf, Mat52, Mat32, Exponential"
  END IF
  
END SUBROUTINE

! NetCDF status check
SUBROUTINE check(stat)

  INTEGER, INTENT(IN) :: stat

  IF(stat /= NF90_NOERR) then
    PRINT *, TRIM(NF90_STRERROR(stat))
    STOP "Stopped"
  END IF

END SUBROUTINE

! Make predictions at new x points
SUBROUTINE predict()
  CONTINUE
END SUBROUTINE

! Kernel functions
FUNCTION rbf()
  REAL(dp) :: rbf
  rbf = 0.0_dp
END FUNCTION
FUNCTION matern52()
  REAL(dp) :: matern52
  matern52 = 0.0_dp
END FUNCTION
FUNCTION matern32()
  REAL(dp) :: matern32
  matern32 = 0.0_dp
END FUNCTION
FUNCTION exponential()
  REAL(dp) :: exponential
  exponential = 0.0_dp
END FUNCTION

! x conversion/reversion functions
! only std_uniform for now

END MODULE

PROGRAM main

  USE faerie

  CALL read_data()

  DEALLOCATE(gp%l,gp%x,gp%y,dist_bnds)

END PROGRAM
