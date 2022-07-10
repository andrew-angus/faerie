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
  FUNCTION kernel(l,var)
    USE ISO_FORTRAN_ENV
    REAL(REAL64), INTENT(IN) :: var
    REAL(REAL64), INTENT(IN), DIMENSION(:) :: l
    REAL(REAL64) :: kernel
  END FUNCTION
END INTERFACE

! GP object
! single y output, zero prior mean, for now
TYPE gp_obj
  REAL(dp) :: noise, var
  REAL(dp), DIMENSION(:), ALLOCATABLE :: l
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: x
  REAL(dp), DIMENSION(:), ALLOCATABLE :: y
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

CONTAINS

! Read in dataset, distributions, and gp params
SUBROUTINE read_data()

  ! Known labels
  CHARACTER (LEN = *), PARAMETER :: fname = "gp.nc"
  !CHARACTER (LEN = *), DIMENSION(3), PARAMETER :: &
  !  attlabs = (/"kernel","var","noise"/)
  !CHARACTER (LEN = *), DIMENSION(4), PARAMETER :: &
  !  dimlabs = (/"inputs","outputs","samples","bounds"/)
  !CHARACTER (LEN = *), DIMENSION(4), PARAMETER :: &
  !  varlabs = (/"lengthscales","input_samples","output_samples", &
  !              "dist_bounds","dists"/)

  ! IDs
  INTEGER :: fid,attlen
  !INTEGER, DIMENSION(3) :: attids
  !INTEGER, DIMENSION(5) :: varids

  ! Outputs
  CHARACTER(:),ALLOCATABLE :: kern

  ! We are reading 2D data, a 6 x 12 grid.
  !integer, parameter :: NX = 6, NY = 12
  !integer :: data_in(NY, NX)

  ! Open file
  CALL check(NF90_OPEN(fname, NF90_NOWRITE, fid))

  ! Get global attributes
  CALL check(NF90_INQUIRE_ATTRIBUTE(fid,NF90_GLOBAL,"kernel",len=attlen))
  ALLOCATE(CHARACTER(attlen)::kern)
  !PRINT *, attlen
  CALL check(NF90_GET_ATT(fid,NF90_GLOBAL,"kernel",kern))

  ! Get variable IDs
  !call check(NF90_INQ_VARID(fid, "data", varid) )

  ! Read the data.
  !call check(NF90_GET_VAR(fid, varid, data_in) )

  ! Close the file, freeing all resources.
  call check(NF90_CLOSE(fid))

  ! Assign GP function pointer
  PRINT *, kern
  
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

END PROGRAM
