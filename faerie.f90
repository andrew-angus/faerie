MODULE faerie

USE ISO_C_BINDING
USE ISO_FORTRAN_ENV

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
  USE ISO_FORTRAN_ENV
  FUNCTION kernel(l,var)
    REAL(REAL64), INTENT(IN) :: var
    REAL(REAL64), INTENT(IN), DIMENSION(:) :: l
    REAL(REAL64) :: kernel
  END FUNCTION
END INTERFACE

! GP object
TYPE gp_obj
  REAL(dp) :: noise, var
  REAL(dp), DIMENSION(:), ALLOCATABLE :: l
  PROCEDURE(kernel), POINTER, NOPASS :: kern => NULL()
END TYPE
TYPE(gp_obj) :: gp

! Read in dataset and hyper-parameters
SUBROUTINE read_data()
  CONTINUE
END SUBROUTINE

! Make predictions at new x points
SUBROUTINE predict()
  CONTINUE
END SUBROUTINE

! Kernel functions
FUNCTION rbf()
  REAL(dp) :: rbf
  CONTINUE
END FUNCTION
FUNCTION matern52()
  REAL(dp) :: matern52
  CONTINUE
END FUNCTION
FUNCTION matern32()
  REAL(dp) :: matern32
  CONTINUE
END FUNCTION
FUNCTION exponential()
  REAL(dp) :: exponential
  CONTINUE
END FUNCTION
