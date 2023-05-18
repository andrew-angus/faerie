MODULE faerie

USE NETCDF
USE IEEE_ARITHMETIC
USE ISO_FORTRAN_ENV

IMPLICIT NONE

SAVE

! Parameters
INTEGER, PARAMETER :: num = REAL64
REAL(num), PARAMETER :: rat13 = 1.0_num/3.0_num
REAL(num), PARAMETER :: sqrt2 = SQRT(2.0_num)
REAL(num), PARAMETER :: sqrt3 = SQRT(3.0_num)
REAL(num), PARAMETER :: sqrt5 = SQRT(5.0_num)
REAL(num), PARAMETER :: pi_num = ACOS(-1.0_num)
REAL(num), PARAMETER :: isqrtpi = 1.0_num/SQRT(pi_num)
! Safety factors on maxmin to allow prediction slightly outside bounds
! and stabilise subsequent kumaraswamy transforms
REAL(num), PARAMETER :: safety = 0.01_num
REAL(num), PARAMETER :: safety2 = 1.0_num/(1.0_num-2.0_num*safety)
! Gauss quadrature points and weights for prediction reversion
REAL(num), DIMENSION(8), PARAMETER :: xgh = &
  (/ -2.930637420257244_num, -1.981656756695843_num, -1.15719371244678_num, &
 -0.381186990207322_num,  0.381186990207322_num,  1.15719371244678_num, &
  1.981656756695843_num,  2.930637420257244_num /)
REAL(num), DIMENSION(8), PARAMETER :: wgh = &
  (/ 1.996040722113678e-04_num, 1.707798300741347e-02_num, &
  2.078023258148918e-1_num, 6.611470125582415e-01_num, &
  6.611470125582415e-01_num, 2.078023258148918e-1_num, &
  1.707798300741347e-02_num, 1.996040722113678e-04_num /)

! Kernel function interface
INTERFACE
  FUNCTION kernel(x1,x2,il,var,nx,beta)
    USE ISO_FORTRAN_ENV
    INTEGER, PARAMETER :: num = REAL64
    REAL(num), DIMENSION(:), INTENT(IN) :: x1,x2,il,var,beta
    INTEGER, INTENT(IN) :: nx
    REAL(num) :: kernel
  END FUNCTION
END INTERFACE

! Conversion/reverson function interfaces
INTERFACE
  FUNCTION convert(a,args)
    USE ISO_FORTRAN_ENV
    INTEGER, PARAMETER :: num = REAL64
    REAL(num), INTENT(IN) :: a
    REAL(num), DIMENSION(:), INTENT(IN) :: args
    REAL(num) :: convert
  END FUNCTION
END INTERFACE
INTERFACE
  FUNCTION revert(a,args)
    USE ISO_FORTRAN_ENV
    INTEGER, PARAMETER :: num = REAL64
    REAL(num), INTENT(IN) :: a
    REAL(num), DIMENSION(:), INTENT(IN) :: args
    REAL(num) :: revert
  END FUNCTION
END INTERFACE

! Mean function interface
INTERFACE
  FUNCTION meanf(x,args)
    USE ISO_FORTRAN_ENV
    INTEGER, PARAMETER :: num = REAL64
    REAL(num), DIMENSION(:), INTENT(IN) :: x, args
    REAL(num) :: meanf
  END FUNCTION
END INTERFACE

! Conversion/reversion class
TYPE conrev
  PROCEDURE(convert), POINTER, NOPASS :: conv => NULL() ! Conversion fun pointer
  PROCEDURE(revert), POINTER, NOPASS :: rev => NULL() ! Reversion fun pointer
END TYPE

! GP object
! single y output, for now
TYPE gp_obj
  TYPE(conrev) :: yconrev ! y conversion/reversions
  REAL(num), DIMENSION(:), ALLOCATABLE :: yconrevargs ! y conversion arguments
  REAL(num), DIMENSION(:), ALLOCATABLE :: y ! Output samples
  REAL(num), DIMENSION(:), ALLOCATABLE :: yc ! Normalised output samples
  PROCEDURE(kernel), POINTER, NOPASS :: kern => NULL() ! Kernel function pointer
  INTEGER :: nsamps ! Number of samples
  TYPE(conrev), DIMENSION(:), ALLOCATABLE :: xconrev ! x conversion/reversions
  REAL(num), DIMENSION(:), ALLOCATABLE :: xconrevargs ! x conversion arguments
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: x ! Input samples
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: xc ! Normalised input samples
  REAL(num), DIMENSION(:), ALLOCATABLE :: il ! Inverse lengthscales
  REAL(num), DIMENSION(:), ALLOCATABLE :: var ! Variances
  REAL(num), DIMENSION(:), ALLOCATABLE :: beta ! Additional kerner parameter
  INTEGER :: nx ! Number of inputs
  REAL(num), DIMENSION(:), ALLOCATABLE :: lowtri, alpha ! Saved intermediates
  REAL(num), DIMENSION(:), ALLOCATABLE :: xmin, xmax ! x data min/max
  PROCEDURE(meanf), POINTER, NOPASS :: mean => NULL() ! mean function pointer
  REAL(num), DIMENSION(:), ALLOCATABLE :: meancoeffs ! coefficients for lr mean
  REAL(num) :: noise ! Gaussian noise, kernel variance
  LOGICAL, DIMENSION(:), ALLOCATABLE :: xnorm ! Whether to meanstd normalise x
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: xmaxmins ! Maximum and minimum samples
  REAL(num), DIMENSION(:), ALLOCATABLE :: maxmins ! Maximum and minimum values
  REAL(num) :: maxvar
  LOGICAL, DIMENSION(:), ALLOCATABLE :: fmaxmins ! Maximum and minimum flags
END TYPE

! Loop counters
INTEGER :: fae,rie

CONTAINS

! Read in dataset, distributions, and gp params from NetCDF file
SUBROUTINE read_gp(fname,gp)

  ! Arguments
  CHARACTER (LEN = *), INTENT(IN) :: fname
  TYPE(gp_obj), INTENT(OUT) :: gp

  ! Known labels
  CHARACTER (LEN = *), DIMENSION(9), PARAMETER :: &
    dimlabs = (/"inputs ","lscales","varns  ","outputs","samples", &
                "maxstr ","xcargs ","ycargs ","mcoeffs"/)
  CHARACTER (LEN = *), DIMENSION(9), PARAMETER :: &
    varlabs = (/"lengthscales  ","variances     ","input_samples ",&
                "output_samples","xconrevs      ","xnorms        ",&
                "xconrevargs   ","yconrevargs   ","meancoeffs    "/)

  ! IDs and lengths
  INTEGER :: fid,attlen
  INTEGER, DIMENSION(9) :: dimids, dimlens, varids

  ! Local output
  INTEGER(1), DIMENSION(:), ALLOCATABLE :: xnorm
  CHARACTER(:),ALLOCATABLE :: kern, yconrev, mean
  REAL(num), DIMENSION(:,:), ALLOCATABLE :: y2d
  CHARACTER(len=30), DIMENSION(:), ALLOCATABLE :: xconrev

  ! Open file
  CALL check(NF90_OPEN(fname, NF90_NOWRITE, fid))

  ! Get global attributes
  CALL check(NF90_INQUIRE_ATTRIBUTE(fid,NF90_GLOBAL,"kernel",len=attlen))
  ALLOCATE(CHARACTER(attlen)::kern)
  CALL check(NF90_INQUIRE_ATTRIBUTE(fid,NF90_GLOBAL,"yconrev",len=attlen))
  ALLOCATE(CHARACTER(attlen)::yconrev)
  CALL check(NF90_INQUIRE_ATTRIBUTE(fid,NF90_GLOBAL,"mean",len=attlen))
  ALLOCATE(CHARACTER(attlen)::mean)
  CALL check(NF90_GET_ATT(fid,NF90_GLOBAL,"kernel",kern))
  CALL check(NF90_GET_ATT(fid,NF90_GLOBAL,"noise",gp%noise))
  CALL check(NF90_GET_ATT(fid,NF90_GLOBAL,"yconrev",yconrev))
  CALL check(NF90_GET_ATT(fid,NF90_GLOBAL,"mean",mean))

  ! Get dimension/variable IDs and lengths
  DO fae = 1, 9
    CALL check(NF90_INQ_DIMID(fid,TRIM(dimlabs(fae)),dimids(fae)))
    CALL check(NF90_INQUIRE_DIMENSION(fid,dimids(fae),len=dimlens(fae)))
    CALL check(NF90_INQ_VARID(fid,TRIM(varlabs(fae)),varids(fae)))
  END DO

  ! Allocate and read variables
  ALLOCATE(gp%il(dimlens(2)),gp%x(dimlens(1),dimlens(5)))
  ALLOCATE(gp%var(dimlens(3)),y2d(dimlens(4),dimlens(5)))
  ALLOCATE(gp%xc(dimlens(1),dimlens(5)),gp%lowtri(dimlens(5)*(dimlens(5)+1)/2))
  ALLOCATE(gp%y(dimlens(5)),gp%yc(dimlens(5)),gp%alpha(dimlens(5)))
  ALLOCATE(gp%xmin(dimlens(1)),gp%xmax(dimlens(1)),xnorm(dimlens(1)))
  ALLOCATE(gp%xnorm(dimlens(1)),gp%xconrev(dimlens(1)),xconrev(dimlens(1)))
  ALLOCATE(gp%xconrevargs(dimlens(7)),gp%yconrevargs(dimlens(8)))
  ALLOCATE(gp%meancoeffs(dimlens(9)))
  CALL check(NF90_GET_VAR(fid, varids(1), gp%il))
  CALL check(NF90_GET_VAR(fid, varids(2), gp%var))
  CALL check(NF90_GET_VAR(fid, varids(3), gp%x))
  CALL check(NF90_GET_VAR(fid, varids(4), y2d))
  CALL check(NF90_GET_VAR(fid, varids(5), xconrev))
  CALL check(NF90_GET_VAR(fid, varids(6), xnorm))
  CALL check(NF90_GET_VAR(fid, varids(7), gp%xconrevargs))
  CALL check(NF90_GET_VAR(fid, varids(8), gp%yconrevargs))
  CALL check(NF90_GET_VAR(fid, varids(9), gp%meancoeffs))

  ! Close file
  CALL check(NF90_CLOSE(fid))

  ! Assign GP object attributes
  gp%nx = dimlens(1)
  gp%nsamps = dimlens(5)
  gp%il = 1/gp%il
  gp%y = y2d(1,:)

  ! Allocate maxmin arrays
  ALLOCATE(gp%maxmins(2*gp%nx),gp%fmaxmins(2*gp%nx),gp%xmaxmins(gp%nx,2*gp%nx))
  DO fae = 1, gp%nx
    gp%maxmins(2*fae-1) = -HUGE(0.0_num)
    gp%maxmins(2*fae) = HUGE(0.0_num)
  END DO
  gp%fmaxmins = .FALSE.
  gp%xmaxmins = 0.0_num
  gp%maxvar = 0.0_num

  ! Kernel function pointer
  IF (kern == "RBF") THEN
    gp%kern => rbf
  ELSE IF (kern == "Matern52") THEN
    gp%kern => matern52
  ELSE IF (kern == "Matern32") THEN
    gp%kern => matern32
  ELSE IF (kern == "Exponential") THEN
    gp%kern => exponential
  ELSE IF (kern == "Exponential+Matern32") THEN
    gp%kern => expmat32
  ELSE IF (kern == "Exponential+Matern52") THEN
    gp%kern => expmat52
  ELSE IF (kern == "Exponential+RBF") THEN
    gp%kern => exprbf
  ELSE IF (kern == "Matern32+RBF") THEN
    gp%kern => mat32rbf
  ELSE
    STOP &
      "Read kernel name invalid:" // &
      " must be one of RBF, Matern52, Matern32, Exponential, or combination." 
  END IF

  ! x conversion/reversion function pointers
  DO fae = 1, gp%nx
    IF (TRIM(xconrev(fae)) == 'logarithm') THEN
      gp%xconrev(fae)%conv => con_log
    ELSE IF (TRIM(xconrev(fae)) == 'identity') THEN
      gp%xconrev(fae)%conv => con_id
    ELSE IF (TRIM(xconrev(fae)) == 'kumaraswamy') THEN
      gp%xconrev(fae)%conv => con_kumaraswamy
    ELSE
      STOP &
        "Read xconrev name invalid:" // &
        " must be one of identity or logarithm." 
    END IF

    ! x normalisation
    IF (xnorm(fae) == 1) THEN
      gp%xnorm(fae) = .TRUE.
    ELSE
      gp%xnorm(fae) = .FALSE.
    END IF

  END DO

  ! y conversion/reversion function pointers
  IF (yconrev == 'logarithm') THEN
    gp%yconrev%conv => con_log
    gp%yconrev%rev => rev_log
  ELSE IF (yconrev == 'identity') THEN
    gp%yconrev%conv => con_id
    gp%yconrev%rev => rev_id
  ELSE IF (yconrev == 'sal') THEN
    gp%yconrev%conv => con_sal
    gp%yconrev%rev => rev_sal
  ELSE IF (yconrev == 'boxcoxsal') THEN
    gp%yconrev%conv => con_boxcoxsal
    gp%yconrev%rev => rev_boxcoxsal
  ELSE
    STOP &
      "Read yconrev name invalid:" // &
      " must be one of identity, logarithm, sal, boxcoxsal." 
  END IF

  ! y mean function pointers
  IF (mean == 'zero') THEN
    gp%mean => zero_mean
  ELSE IF (mean == 'linear') THEN
    gp%mean => lr_mean
  ELSE
    STOP &
      "Read mean name invalid:" // &
      " must be one of zero or linear." 
  END IF
  
  ! Convert datasets
  CALL full_convert(gp)

  ! Deallocate temp arrays
  DEALLOCATE(kern,yconrev,y2d,xnorm,xconrev)
  
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
SUBROUTINE full_convert(gp)

  TYPE(gp_obj), INTENT(INOUT) :: gp
  REAL(num) :: xmax, xmin, xdiff

  ! x normalisation
  DO fae = 1, gp%nx
    IF (gp%xnorm(fae)) THEN
      xmin = MINVAL(gp%x(fae,:))
      xmax = MAXVAL(gp%x(fae,:))
      xdiff = safety2*(xmax-xmin)
      gp%xmin(fae) = xmin - safety*xdiff
      gp%xmax(fae) = xmax + safety*xdiff
      DO rie = 1, gp%nsamps
        gp%xc(fae,rie) = &
          convert_minmax(gp%x(fae,rie),gp%xmin(fae),gp%xmax(fae))
      END DO
    ELSE
      gp%xc(fae,:) = gp%x(fae,:)
    END IF
  END DO

  ! x conversions
  DO fae = 1, gp%nsamps
    DO rie = 1, gp%nx
      gp%xc(rie,fae) = gp%xconrev(rie)%conv(gp%xc(rie,fae),&
        gp%xconrevargs(2*(rie-1)+1:2*rie))
    END DO
  END DO

  ! y conversion
  DO fae = 1, gp%nsamps
    gp%yc(fae) = gp%yconrev%conv(gp%y(fae) &
      - gp%mean(gp%x(:,fae),gp%meancoeffs) &
      ,gp%yconrevargs)
  END DO

END SUBROUTINE

! Evaluate covariance matrix for input samples dataset
SUBROUTINE data_covariances(gp)

  TYPE(gp_obj), INTENT(INOUT) :: gp
  INTEGER :: cnt

  ! Evaluate kernel at each unique combination of input data vectors
  ! In column order by lower triangle
  cnt = 1
  DO fae = 1, gp%nsamps
    DO rie = fae, gp%nsamps
      IF (fae == rie) THEN
        ! Add noise on diagonal
        gp%lowtri(cnt) = SUM(gp%var) + gp%noise
      ELSE
        gp%lowtri(cnt) = &
          gp%kern(gp%xc(:,fae),gp%xc(:,rie),gp%il,gp%var,gp%nx,gp%beta) 
      END IF
      cnt = cnt + 1
    END DO
  END DO

END SUBROUTINE

! Use LAPACK Cholesky decomposition solver to obtain
! alpha = (K+v_nI)^-1.y and L (where K+v_nI = LL^T)
! For a fixed dataset these are also fixed
SUBROUTINE cholesky_solve(gp)

  TYPE(gp_obj), INTENT(INOUT) :: gp
  INTEGER :: info

  ! Copy alpha array with y to be overwritten 
  gp%alpha = gp%yc

  ! Obtain lower triangle of symmetric data covariance matrix
  CALL data_covariances(gp)
  
  ! LAPACK Cholesky solve
  CALL dppsv('L',gp%nsamps,1,gp%lowtri,gp%alpha,gp%nsamps,info)

END SUBROUTINE

! Make predictions at new x points
SUBROUTINE predict(gp,x,ypmean,ypvar)

  TYPE(gp_obj), INTENT(IN) :: gp
  !REAL(num), DIMENSION(:), INTENT(IN) :: x
  REAL(num), DIMENSION(:,:), INTENT(IN) :: x
  REAL(num), DIMENSION(:), INTENT(OUT) :: ypmean, ypvar

  REAL(num), DIMENSION(gp%nx) :: xconv
  REAL(num), DIMENSION(gp%nsamps) :: kstar
  REAL(num), DIMENSION(8) :: ygh
  INTEGER :: sx

  ! Loop over prediction points
  sx = SIZE(x,2)
  DO fae = 1, sx

    ! Convert x
    DO rie = 1, gp%nx
      IF (gp%xnorm(rie)) THEN
        xconv(rie) = convert_minmax(x(rie,fae),gp%xmin(rie),gp%xmax(rie))
      END IF
      xconv(rie) = gp%xconrev(rie)%conv(xconv(rie),&
        gp%xconrevargs(2*rie-1:2*rie))
    END DO

    ! Get kstar
    DO rie = 1, gp%nsamps
      kstar(rie) = gp%kern(xconv,gp%xc(:,rie),gp%il,gp%var,gp%nx,gp%beta)
    END DO

    ! Get mean and variance prediction
    ypmean(fae) = DOT_PRODUCT(kstar,gp%alpha)

    ! LAPACK solve for variance vector
    CALL dtpsv('L','N','N',gp%nsamps,gp%lowtri,kstar,1)

    ! Get variance prediction
    ypvar(fae) = SUM(gp%var) - DOT_PRODUCT(kstar,kstar)

    ! Revert y predictions
    ! Gauss quadrature, see Rios, G. & Tobar, F.
    ! Compositionally-warped Gaussian processes. 
    ! Neural Networks 118, 235–246 (2019).
    ygh = sqrt2*SQRT(ypvar(fae))*xgh+ypmean(fae)
    DO rie = 1, 8
      ygh(rie) = gp%yconrev%rev(ygh(rie),gp%yconrevargs)
    END DO
    ygh = ygh + gp%mean(x(:,fae),gp%meancoeffs)
    ypmean(fae) = isqrtpi*SUM(wgh*ygh)
    ypvar(fae) = isqrtpi*SUM(wgh*ygh**2) - ypmean(fae)**2

    ! Check on variance prediction
    IF ((IEEE_IS_NAN(ypvar(fae))) .OR. (ypvar(fae) > HUGE(ypvar(fae)))) THEN
      PRINT *, "Error: Positive definite matrix not" // &
        " produced by kernel"
      STOP
    END IF

  END DO

END SUBROUTINE

! Kernel functions
! Radial basis function kernel
FUNCTION rbf(x1,x2,il,var,nx,beta)
  REAL(num), DIMENSION(:), INTENT(IN) :: x1,x2,il,var,beta
  INTEGER, INTENT(IN) :: nx
  REAL(num), DIMENSION(nx) :: r
  REAL(num) :: rbf
  r = (x1 - x2)*il
  rbf = var(1)*EXP(-0.5_num*DOT_PRODUCT(r,r))
END FUNCTION
! Matern52
FUNCTION matern52(x1,x2,il,var,nx,beta)
  REAL(num), DIMENSION(:), INTENT(IN) :: x1,x2,il,var,beta
  INTEGER, INTENT(IN) :: nx
  REAL(num), DIMENSION(nx) :: r
  REAL(num) :: matern52, ril
  r = (x1 - x2)*il
  ril = sqrt5*SQRT(DOT_PRODUCT(r,r))
  matern52 = var(1)*(1+ril+rat13*ril**2)*EXP(-ril)
END FUNCTION
! Matern32
FUNCTION matern32(x1,x2,il,var,nx,beta)
  REAL(num), DIMENSION(:), INTENT(IN) :: x1,x2,il,var,beta
  INTEGER, INTENT(IN) :: nx
  REAL(num), DIMENSION(nx) :: r
  REAL(num) :: matern32, ril
  r = (x1 - x2)*il
  ril = sqrt3*SQRT(DOT_PRODUCT(r,r))
  matern32 = var(1)*(1+ril)*EXP(-ril)
END FUNCTION
! Exponential
FUNCTION exponential(x1,x2,il,var,nx,beta)
  REAL(num), DIMENSION(:), INTENT(IN) :: x1,x2,il,var,beta
  INTEGER, INTENT(IN) :: nx
  REAL(num), DIMENSION(nx) :: r
  REAL(num) :: exponential
  r = (x1 - x2)*il
  exponential = var(1)*EXP(-SQRT(DOT_PRODUCT(r,r)))
END FUNCTION

! Kernel combinations
! exponential + matern32
FUNCTION expmat32(x1,x2,il,var,nx,beta)
  REAL(num), DIMENSION(:), INTENT(IN) :: x1,x2,il,var,beta
  INTEGER, INTENT(IN) :: nx
  REAL(num) :: expmat32
  expmat32 = exponential(x1,x2,il(1:nx),(/var(1)/),nx,beta) &
    + matern32(x1,x2,il(nx+1:nx*2),(/var(2)/),nx,beta)
END FUNCTION
! exponential + matern52
FUNCTION expmat52(x1,x2,il,var,nx,beta)
  REAL(num), DIMENSION(:), INTENT(IN) :: x1,x2,il,var,beta
  INTEGER, INTENT(IN) :: nx
  REAL(num) :: expmat52
  expmat52 = exponential(x1,x2,il(1:nx),(/var(1)/),nx,beta) &
    + matern52(x1,x2,il(nx+1:nx*2),(/var(2)/),nx,beta)
END FUNCTION
! exponential + RBF
FUNCTION exprbf(x1,x2,il,var,nx,beta)
  REAL(num), DIMENSION(:), INTENT(IN) :: x1,x2,il,var,beta
  INTEGER, INTENT(IN) :: nx
  REAL(num) :: exprbf
  exprbf = exponential(x1,x2,il(1:nx),(/var(1)/),nx,beta) &
    + rbf(x1,x2,il(nx+1:nx*2),(/var(2)/),nx,beta)
END FUNCTION
! matern32 + RBF
FUNCTION mat32rbf(x1,x2,il,var,nx,beta)
  REAL(num), DIMENSION(:), INTENT(IN) :: x1,x2,il,var,beta
  INTEGER, INTENT(IN) :: nx
  REAL(num) :: mat32rbf
  mat32rbf = matern32(x1,x2,il(1:nx),(/var(1)/),nx,beta) &
    + rbf(x1,x2,il(nx+1:nx*2),(/var(2)/),nx,beta)
END FUNCTION

! x Normalisation by data min and max
FUNCTION convert_minmax(a,amin,amax)
  REAL(num), INTENT(IN) :: a, amin, amax
  REAL(num) :: convert_minmax
  convert_minmax = (a-amin)/(amax-amin)
END FUNCTION

! x/y conversion/reversion routines
! identity
FUNCTION con_id(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: con_id
  con_id = a
END FUNCTION
FUNCTION rev_id(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: rev_id
  rev_id = a
END FUNCTION
! affine
FUNCTION con_affine(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: con_affine
  con_affine = args(1) + args(2)*a
END FUNCTION
FUNCTION rev_affine(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: rev_affine
  rev_affine = (a - args(1))/args(2)
END FUNCTION
! kumaraswamy (conversion only)
FUNCTION con_kumaraswamy(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num) :: con_kumaraswamy
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  con_kumaraswamy = 1.0_num - (1.0_num - a**args(1))**args(2)
END FUNCTION
! logarithm
FUNCTION con_log(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: con_log
  con_log = LOG(a)
END FUNCTION
FUNCTION rev_log(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: rev_log
  rev_log = EXP(a)
END FUNCTION
! sinharcsinh
FUNCTION con_sarc(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: con_sarc
  con_sarc = SINH(args(2)*ASINH(a)-args(1))
END FUNCTION
FUNCTION rev_sarc(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: rev_sarc
  rev_sarc = SINH((ASINH(a)+args(1))/args(2))
END FUNCTION
! boxcox
FUNCTION con_boxcox(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: con_boxcox
  con_boxcox = (SIGN(1.0_num,a)*ABS(a)**args(1)-1)/args(1)
END FUNCTION
FUNCTION rev_boxcox(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: rev_boxcox, yterm
  yterm = a*args(1) + 1
  rev_boxcox = SIGN(1.0_num,yterm)*ABS(yterm)**(1/args(1))
END FUNCTION

! Standard multi-layer transforms
FUNCTION con_sal(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: con_sal
  con_sal = con_affine(a,args(1:2))
  con_sal = con_sarc(con_sal,args(3:4))
  con_sal = con_affine(con_sal,args(5:6))
  con_sal = con_sarc(con_sal,args(7:8))
  con_sal = con_affine(con_sal,args(9:10))
END FUNCTION
FUNCTION rev_sal(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: rev_sal
  rev_sal = rev_affine(a,args(9:10))
  rev_sal = rev_sarc(rev_sal,args(7:8))
  rev_sal = rev_affine(rev_sal,args(5:6))
  rev_sal = rev_sarc(rev_sal,args(3:4))
  rev_sal = rev_affine(rev_sal,args(1:2))
END FUNCTION
FUNCTION con_boxcoxsal(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: con_boxcoxsal
  con_boxcoxsal = con_boxcox(a,args(1:1))
  con_boxcoxsal = con_sal(con_boxcoxsal,args(2:11))
END FUNCTION
FUNCTION rev_boxcoxsal(a,args)
  REAL(num), INTENT(IN) :: a
  REAL(num), DIMENSION(:), INTENT(IN) :: args
  REAL(num) :: rev_boxcoxsal
  rev_boxcoxsal = rev_sal(a,args(2:11))
  rev_boxcoxsal = rev_boxcox(rev_boxcoxsal,args(1:1))
END FUNCTION

! Mean functions
! zero
FUNCTION zero_mean(x,args)
  REAL(num), DIMENSION(:), INTENT(IN) :: x, args
  REAL(num) :: zero_mean
  zero_mean = 0.0_num
END FUNCTION
! linear regression
FUNCTION lr_mean(x,args)
  REAL(num), DIMENSION(:), INTENT(IN) :: x, args
  REAL(num) :: lr_mean
  lr_mean = args(1) + DOT_PRODUCT(x,args(2:))
END FUNCTION

END MODULE
