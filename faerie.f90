MODULE faerie

USE ISO_C_BINDING
USE ISO_FORTRAN_ENV
USE NETCDF

IMPLICIT NONE

SAVE

! Parameters
INTEGER, PARAMETER :: cd = C_DOUBLE
INTEGER, PARAMETER :: ci = C_INT
INTEGER, PARAMETER :: cb = C_BOOL
INTEGER, PARAMETER :: dp = REAL64
INTEGER, PARAMETER :: fi = INT32
REAL(dp), PARAMETER :: pi = ACOS(-1.0_dp)
REAL(dp), PARAMETER :: small = 1e-30_dp
REAL(dp), PARAMETER :: releps = 1.0e-6_dp
REAL(dp), PARAMETER :: rat53 = 5.0_dp/3.0_dp
REAL(dp), PARAMETER :: rat13 = 1.0_dp/3.0_dp
REAL(dp), PARAMETER :: sqrt3 = SQRT(3.0_dp)
REAL(dp), PARAMETER :: sqrt5 = SQRT(5.0_dp)
CHARACTER(LEN=1), PARAMETER :: creturn = ACHAR(13)

! Kernel function interface
INTERFACE
  FUNCTION kernel(x1,x2)
    USE ISO_FORTRAN_ENV
    REAL(REAL64), DIMENSION(:), INTENT(IN) :: x1,x2
    REAL(REAL64) :: kernel
  END FUNCTION
END INTERFACE

! GP object
! single y output, zero prior mean, for now
TYPE gp_obj
  INTEGER(fi) :: nx,nsamps ! Number of inputs and samples
  PROCEDURE(kernel), POINTER, NOPASS :: kern => NULL() ! Kernel function pointer
  REAL(dp) :: noise, var ! Gaussian noise, kernel variance
  REAL(dp) :: ymean, ystd ! Mean and std dev of y data for normalisation
  REAL(dp), DIMENSION(:), ALLOCATABLE :: il ! Inverse lengthscales
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: x ! Input samples
  REAL(dp), DIMENSION(:), ALLOCATABLE :: y ! Output samples
  REAL(dp), DIMENSION(:,:), ALLOCATABLE :: xc ! Normalised input samples
  REAL(dp), DIMENSION(:), ALLOCATABLE :: yc ! Normalised output samples
  REAL(dp), DIMENSION(:), ALLOCATABLE :: xmean, xstd ! x data mean/ std dev
END TYPE
TYPE(gp_obj) :: gp

! Shared data
INTEGER(fi) :: i,j,info
REAL(dp), DIMENSION(:,:), ALLOCATABLE :: y2d
REAL(dp), DIMENSION(:), ALLOCATABLE :: lowtri, alpha
REAL(dp) :: ymean, yvar

CONTAINS

! Read in dataset, distributions, and gp params from NetCDF file
SUBROUTINE read_data(fname)

  ! Known labels
  CHARACTER (LEN = *), INTENT(IN) :: fname
  CHARACTER (LEN = *), DIMENSION(3), PARAMETER :: &
    dimlabs = (/"inputs ","outputs","samples"/)
  CHARACTER (LEN = *), DIMENSION(3), PARAMETER :: &
    varlabs = (/"lengthscales  ","input_samples ","output_samples"/)

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
  DO i = 1, 3
    CALL check(NF90_INQ_DIMID(fid,TRIM(dimlabs(i)),dimids(i)))
    CALL check(NF90_INQUIRE_DIMENSION(fid,dimids(i),len=dimlens(i)))
    CALL check(NF90_INQ_VARID(fid,TRIM(varlabs(i)),varids(i)))
  END DO

  ! Allocate and read variables
  ALLOCATE(gp%il(dimlens(1)),gp%x(dimlens(3),dimlens(1)),y2d(dimlens(3),dimlens(2)))
  ALLOCATE(gp%xc(dimlens(3),dimlens(1)),lowtri(dimlens(3)*(dimlens(3)+1)/2))
  ALLOCATE(gp%y(dimlens(3)),gp%yc(dimlens(3)),alpha(dimlens(3)))
  ALLOCATE(gp%xmean(dimlens(1)),gp%xstd(dimlens(1)))
  CALL check(NF90_GET_VAR(fid, varids(1), gp%il))
  CALL check(NF90_GET_VAR(fid, varids(2), gp%x))
  CALL check(NF90_GET_VAR(fid, varids(3), y2d))

  ! Close file
  CALL check(NF90_CLOSE(fid))

  ! Assign GP function pointer and dimensions
  gp%nx = dimlens(1)
  gp%nsamps = dimlens(3)
  gp%il = 1/gp%il
  gp%y = y2d(:,1)
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

  ! Convert datasets
  CALL full_convert()

  ! Deallocate temp arrays
  DEALLOCATE(kern,y2d)
  
END SUBROUTINE

! NetCDF status check
SUBROUTINE check(stat)

  INTEGER, INTENT(IN) :: stat

  IF(stat /= NF90_NOERR) then
    PRINT *, TRIM(NF90_STRERROR(stat))
    STOP "Stopped"
  END IF

END SUBROUTINE

! Convert full datasets and assign to GP object
SUBROUTINE full_convert()

  ! Get means and standard deviations
  gp%ymean = mean(gp%y)
  gp%ystd = std_dev(gp%y,gp%ymean)
  DO i = 1, gp%nx
    gp%xmean(i) = mean(gp%x(:,i))
    gp%xstd(i) = std_dev(gp%x(:,i),gp%xmean(i))
  END DO

  ! y conversions
  DO j = 1, gp%nsamps
    gp%yc(j) = convert(gp%y(j),gp%ymean,gp%ystd)
  END DO

  ! x-conversions
  DO i = 1, gp%nx
    DO j = 1, gp%nsamps
      gp%xc(j,i) = convert(gp%x(j,i),gp%xmean(i),gp%xstd(i))
    END DO
  END DO

END SUBROUTINE

! Evaluate covariance matrix for input samples dataset
SUBROUTINE data_covariances()

  INTEGER(fi) :: cnt

  ! Evaluate kernel at each unique combination of input data vectors
  ! In column order by lower triangle
  cnt = 1
  DO i = 1, gp%nsamps
    DO j = i, gp%nsamps
      lowtri(cnt) = gp%kern(gp%xc(i,:),gp%xc(j,:)) 
      IF (i == j) THEN
        ! Add noise on diagonal
        lowtri(cnt) = lowtri(cnt) + gp%noise
      END IF
      cnt = cnt + 1
    END DO
  END DO

END SUBROUTINE

! Use LAPACK Cholesky decomposition solver to obtain
! alpha = (K+v_nI)^-1.y and L (where K+v_nI = LL^T)
! For a fixed dataset these are also fixed
SUBROUTINE cholesky_solve()

  ! Copy alpha array with y to be overwritten 
  alpha = gp%yc

  ! Obtain lower triangle of symmetric data covariance matrix
  CALL data_covariances()
  
  ! LAPACK Cholesky solve
  CALL dppsv('L',gp%nsamps,1,lowtri,alpha,gp%nsamps,info)

END SUBROUTINE

! Make predictions at new x points
! Single x vector for now
SUBROUTINE predict(x,ypmean,ypvar)

  REAL(dp), DIMENSION(:), INTENT(INOUT) :: x
  REAL(dp), INTENT(OUT) :: ypmean, ypvar

  REAL(dp), DIMENSION(gp%nsamps) :: kstar

  ! Convert x
  DO i = 1, gp%nx
    x(i) = convert(x(i),gp%xmean(i),gp%xstd(i))
  END DO

  ! Get kstar
  DO i = 1, gp%nsamps
    kstar(i) = gp%kern(gp%xc(i,:),x)
  END DO

  ! Get mean prediction
  ypmean = DOT_PRODUCT(kstar,alpha)

  ! LAPACK solve for variance vector
  CALL dtpsv('L','N','N',gp%nsamps,lowtri,kstar,1)
  
  ! Get variance prediction
  ypvar = gp%kern(x,x) - DOT_PRODUCT(kstar,kstar)

  ! Revert y predictions
  ypmean = revert(ypmean,gp%ymean,gp%ystd)
  ypvar = ypvar*gp%ystd**2

END SUBROUTINE

! Kernel functions
FUNCTION rbf(x1,x2)
  REAL(dp), DIMENSION(:), INTENT(IN) :: x1,x2
  REAL(dp), DIMENSION(gp%nx) :: r
  REAL(dp) :: rbf
  r = (x1 - x2)*gp%il
  rbf = gp%var*EXP(-0.5_dp*DOT_PRODUCT(r,r))
END FUNCTION
FUNCTION matern52(x1,x2)
  REAL(dp), DIMENSION(:), INTENT(IN) :: x1,x2
  REAL(dp), DIMENSION(gp%nx) :: r
  REAL(dp) :: matern52, ril
  r = (x1 - x2)*gp%il
  ril = sqrt5*SQRT(DOT_PRODUCT(r,r))
  matern52 = gp%var*(1+ril+rat13*ril**2)*EXP(-ril)
END FUNCTION
FUNCTION matern32(x1,x2)
  REAL(dp), DIMENSION(:), INTENT(IN) :: x1,x2
  REAL(dp), DIMENSION(gp%nx) :: r
  REAL(dp) :: matern32, ril
  r = (x1 - x2)*gp%il
  ril = sqrt5*SQRT(DOT_PRODUCT(r,r))
  matern32 = gp%var*(1+ril)*EXP(-ril)
END FUNCTION
FUNCTION exponential(x1,x2)
  REAL(dp), DIMENSION(:), INTENT(IN) :: x1,x2
  REAL(dp), DIMENSION(gp%nx) :: r
  REAL(dp) :: exponential, ril
  r = (x1 - x2)*gp%il
  ril = sqrt5*SQRT(DOT_PRODUCT(r,r))
  exponential = gp%var*EXP(-ril)
END FUNCTION

! x/y conversion/reversion routines
! only normalisation by mean and std deviation for now
FUNCTION convert(a,amean,astd)

  REAL(dp), INTENT(IN) :: a, amean, astd
  REAL(dp) :: convert

  convert = (a - amean)/astd

END FUNCTION
FUNCTION revert(a,amean,astd)

  REAL(dp), INTENT(IN) :: a, amean, astd
  REAL(dp) :: revert

  revert = a*astd + amean

END FUNCTION
! Mean and std dev functions
FUNCTION mean(a)

  REAL(dp), DIMENSION(:), INTENT(IN) :: a
  REAL(dp) :: mean

  mean = SUM(a)/SIZE(a)

END FUNCTION
FUNCTION std_dev(a,amean)

  REAL(dp), DIMENSION(:), INTENT(IN) :: a
  REAL(dp), INTENT(IN) :: amean
  REAL(dp) :: std_dev

  std_dev = SQRT(SUM((a-amean)**2)/SIZE(a))

END FUNCTION


END MODULE

! Test program
PROGRAM main

  USE faerie

  REAL(dp), DIMENSION(3) :: &
    xnew = (/7.63838986e+07,1.90197926e+27,3.28797069e+15/)
    !xnew = (/0.48058382,1.01261404,-1.00189514/)

  CALL read_data('gp.nc')

  CALL cholesky_solve()

  CALL predict(xnew,ymean,yvar)
  PRINT *, ymean, yvar

  DEALLOCATE(gp%il,gp%x,gp%y,gp%xc,gp%yc,lowtri,alpha)
  DEALLOCATE(gp%xmean,gp%xstd)

END PROGRAM
