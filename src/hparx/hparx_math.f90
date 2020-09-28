!-*-f90-*-
 
! @LICENSE_HEADER_START@
!
!   This file is part of HPARX.
!   
!   --
!   HPARX: Fortran code library for atmospheric sciences
!   
!   Copyright (C) 2006-2015
!   Hironobu Iwabuchi, Souichiro Hioki, and Rintaro Okamura
!   
!   HPARX is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!   
!   HPARX is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with HPARX. If not, see <http://www.gnu.org/licenses/>.
!
! @LICENSE_HEADER_END@
 




!+
! Library of mathematical utilities
!-
module hparx_math

  use globals
  use hparx_lina, only : la_Gauss_solve
  use hparx_vecmat, only : gridIdx_loc
  implicit none
  private

  ! Public
  public :: cycleRange_R
  public :: effAve_exp_R
  public :: erf_R             ! Error function
  public :: erfc_R            ! Complementary error function
  public :: erfc_cheb_R       ! Complementary error function using a Chebyshev approximation
  public :: fdif_cent_1R      ! finite differences by the central differences
  public :: fitpoly
  public :: gammaLog_R
  public :: gauss_R
  public :: gaussLegen        ! coefficients of Gauss-Legendre quadrature
  public :: grid_refine_1R
  public :: integ_trapz_R
  public :: integ_trapz_eqX_R
  public :: intp_R
  public :: intp_Akima_coef_1R
  public :: intp_bilinear_R     ! bilinear interpolation
  public :: intp_bicubic_eqXY_R ! bicubic interpolation
  public :: intp_cubic_R        ! cubic interpolation
  public :: intp_cubic_eqX_R    ! cubic interpolation
  public :: intp_Hermite3_R
  public :: intp_Hermite3_1R
  public :: intp_Hermite3_coef
  public :: intp_spline_R        ! interpolation of a scalar by natural cubic spline 
  public :: intp_spline_1R       ! interpolation of a vector by natural cubic spline 
  public :: intp_spline_coef_1R  ! 2nd derivatives of a cubic spline
  public :: intp_tabLin_eqX_R
  public :: intp_tabCub_eqX_R
  public :: intpTab_Akima_R
  public :: intpTab_Akima_1R
  public :: intpTab_Akima_oR
  public :: intpTab_Akima_o1R
  public :: intpTab_Lagr_R
  public :: intpTab_linear_R
  public :: intpTab_linear_1R
  public :: intpTab_spline_1R
  public :: legenPoly_1R
  public :: legenPolySer_1R
  public :: modifyTotal_1R
  public :: relaDevMod_1R     ! relative deviation vector modified by a forcing
  public :: root_quad         ! solve a quadratic equation
  public :: root_quadA1       ! solve a quadratic equation with A=1
  public :: root_quadA1_one_R ! a root of quadratic equation
  public :: secAve_exp_1R
  public :: secAve_lin_1R
  public :: smooth_triang_1R
  public :: smooth_noHotSpots ! smooth hot spots in a vector
  public :: sphBesJ_1RD
  public :: sphBesY_1RD

  ! Code shopping list
  !  - special functions (e.g., imcomplete Gamma, Bessel)

contains

  !+
  ! Returns a value reset in a cyclic range
  !-
  function cycleRange_R(x, xlo, xhi) result (res)

    real(R_),  intent(in) :: xlo ! lower bound
    real(R_),  intent(in) :: xhi ! upper bound
    real(R_),  intent(in) :: x   ! a value
    real(R_) :: res ! result

    if (x >= xlo .and. x < xhi) then
       res = x
    else
       res = mod(x - xlo, xhi - xlo)
       if (res < 0.0_R_) then
          res = res + xhi
       else
          res = res + xlo
       end if
    end if

  end function cycleRange_R


  !+
  ! Effective layer average of property with exponential profile exp(beta*x)
  !-
  function effAve_exp_R(b, t) result(a)

    real(R_), intent(in) :: b, t  ! bottom & top values (should be > 0)
    real(R_) :: a, rat

    rat = t / b
    if (b * t <= 0.0_R_ .or. abs(rat - 1.0_R_) < 0.001_R_) then
       a = 0.5_R_ * (b + t)
    else
       a = (t - b) / log(rat)
       !// 03/06/2009, debug, removed a factor (t / b)
    end if

  end function effAve_exp_R


  !+
  ! Error function
  !-
  function erf_R(x) result(y)

    real(R_), intent(in)  :: x    ! X
    real(R_)  :: y    ! Y = erf(X)
    !y = erf(x) ! erf() is nonstandard!
    if (x > 0.0_R_) then
       y = 1.0_R_ - erfc_Cheb_R(x)
    else
       y = erfc_Cheb_R(-x) - 1.0_R_
    end if

  end function erf_R


  !+
  ! Complementary error function
  !-
  function erfc_R(x) result(y)

    real(R_), intent(in)  :: x    ! X
    real(R_)  :: y    ! Y = erfc(X)
    !y = erfc(x) ! erfc() is nonstandard!
    if (x > 0.0_R_) then
       y = erfc_Cheb_R(x)
    else
       y = 2.0_R_ - erfc_Cheb_R(-x)
    end if

  end function erfc_R


  !+
  ! Complementary error function using a Chebyshev approximation (single precision version)
  !  -with relative error less than 1.2e-7
  !-
  function erfc_Cheb_R(z) result(y)

    real(R_), intent(in)  :: z ! abs(X) > 0
    real(R_)  :: y ! Y = erfc(z)
    real(R_)  :: t

    if (z < 0.0_R_) stop 'erfc_Cheb_R: Z should be >= 0.'
    t = 2.0_R_ / (2.0_R_ + z)
    y = -z**2 - 1.26551223_R_ + t*(1.00002368_R_ +t*(0.37409196_R_ + t*(0.09678418_R_ &
         + t*(-0.18628806_R_ + t*(0.27886807_R_ + t*(-1.13520398_R_ + t*(1.4851587_R_ &
         + t*(-0.82215223_R_ + t*0.17087277_R_))))))))
    y = t * exp(y)

  end function erfc_Cheb_R


  !+
  ! Finite difference by central differences
  !-
  function fdif_cent_1R(y, x, dx, mend) result(res)

    real(R_), intent(in) :: y(:) ! y values
    real(R_), intent(in), optional :: x(:) ! x values in a general case
    real(R_), intent(in), optional :: dx   ! x spacing if equally-spaced grid (default = 1.0)
    !// If both of x() and dx are not given, dx = 1.0 is assumed.
    integer,  intent(in), optional :: mend ! method flag for the end points (default = 0)
    !// = -1 : cyclic boundary condition is assumed
    !      0 : 1st-order difference is used
    !      1 : 2nd derivatives at the end points are assumed to be 0
    real(R_) :: res(size(y)) ! result, derivatives
    integer :: m, n, i
    real(R_) :: dold, dnew, d

    ! For the majority
    n = size(y)
    d = 1.0_R_
    if (present(dx)) d = dx
    if (present(x)) then
       dold = (y(2) - y(1)) / (x(2) - x(1))
       do i = 2, n - 1
          dnew = (y(i + 1) - y(i)) / (x(i + 1) - x(i))
          res(i) = 0.5_R_ * (dnew + dold)
          dold = dnew
       end do
    else
       res(2:n-1) = 0.5_R_ / d * (y(3:n) - y(1:n-2))
    end if

    ! End points
    m = 0
    if (present(mend)) m = mend
    if (m == -1) then ! cyclic
       if (present(x)) then
          dold = (y(n) - y(1)) / (x(n) - x(1))
          dnew = (y(2) - y(1)) / (x(2) - x(1))
          res(1) = 0.5_R_ * (dnew + dold)
          dnew = dold
          dold = (y(n) - y(n - 1)) / (x(n) - x(n - 1))
          res(n) = 0.5_R_ * (dnew + dold)
       else
          res(1) = 0.5_R_ / d * (y(2) - y(n))
          res(n) = 0.5_R_ / d * (y(1) - y(n - 1))
       end if
    else ! m = 0 or 1
       if (present(x)) then
          res(1) = (y(2) - y(1)) / (x(2) - x(1)) ! 1st-order difference
          res(n) = (y(n) - y(n - 1)) / (x(n) - x(n - 1))
       else
          res(1) = (y(2) - y(1)) / d
          res(n) = (y(n) - y(n - 1)) / d
       end if
       if (m == 1) then ! 2nd derivatives of zero
          res(1) = 1.5_R_ * res(1) - 0.5_R_ * res(2)
          res(n) = 1.5_R_ * res(n) - 0.5_R_ * res(n - 1)
       end if
    end if

  end function fdif_cent_1R


  !+
  ! Fit a polynomial, Y = a0 + a1*X + a2*X^2 + ...
  !-
  subroutine fitpoly(x, y, coef, rmsey, ierr)

    real(R_), intent(in)  :: x(:), y(:) ! X and Y data
    real(R_), intent(out) :: coef(:)    ! polynomial coefficients
    real(R_), intent(out) :: rmsey      ! RMSE of Y
    integer,  intent(out) :: ierr       ! error code (0=normal, otherwise error)
    real(R_) :: a(size(coef), size(coef)), b(size(coef), 1), fx(size(x))
    integer :: i, j, order

    ! Calculate the coefficients
    order = size(coef) - 1
    do i = 1, order + 1
       do j = 1, order + 1
          a(i, j) = sum(x(:)**(i+j-2))
       end do
       b(i,1) = sum((x(:)**(i-1))*y(:))
    end do
    call la_Gauss_solve(a, b, ierr)
    if (ierr /= 0) return
    coef(:) = b(:,1)

    ! Compute root-mean-square of y residuals
    fx(:) = 0.0_R_
    do i = order + 1, 1, -1
       fx(:) = fx(:) + b(i, 1)*(x(:)**(i-1))
    end do
    rmsey = sqrt(sum((y(:) - fx(:))**2) / size(x))

  end subroutine fitpoly


  !+
  ! Log(Gamma function)
  !-
  function gammaLog_R(x) result(res)

    real(R_),  intent(in) :: x  ! argument (>0) of Gamma function
    real(R_)  :: res ! result
    integer   :: i
    real(RD_) :: sum, tmp, xx
    real(RD_), save :: c0 = 1.000000000190015_RD_
    real(RD_), save :: c(6) = &
         (/ 76.18009172947146_RD_, -86.50532032941677_RD_, 24.01409824083091_RD_, &
         -1.231739572450155_RD_, 0.1208650973866179e-2_RD_, -0.5395239384953e-5_RD_ /)
    real(RD_), save :: fac = 2.5066282746310005_RD_ ! sqrt(2*pi)

    xx = real(x, RD_)
    sum = c0
    do i = 1, 6
       sum = sum + c(i) / (xx + real(i, RD_))
    end do
    tmp = xx + 5.5_RD_
    res = (xx + 0.5_RD_) * log(tmp) - tmp + log(fac * sum / xx) 

  end function gammaLog_R


  !+
  ! Gauss function
  !-
  function gauss_R(x) result(res) 

    real(R_), intent(in) :: x
    real(R_) :: res
    res = exp(-0.5_R_ * x * x) / sqrt(PI2_)

  end function gauss_R


  !+
  ! Get quadrature points and weights of Gauss-Legendre integration
  !  formula for the interval [-1, +1]
  !-
  subroutine gaussLegen(nx, x, w)

    integer,  intent(in)  :: nx   ! # of points
    real(R_), intent(out) :: x(:) ! abcissa in [-1, 1], increasing order
    real(R_), intent(out) :: w(:) ! weights, sum(w(1:nx))=1
    integer  :: nh, ih, i
    real(RD_) :: p0, p1, p2, pp, z, dz

    nh = (nx + 1) / 2
    do ih = 1, nh
       z = cos(DPI_ * (ih - 0.25_RD_) / (nx + 0.5_RD_))
       do       ! Newtonian iteration
          p1 = 0.0_RD_
          p0 = 1.0_RD_
          do i = 1, nx
             p2 = p1
             p1 = p0
             p0 = ((2*i - 1) * z * p1 - (i - 1) * p2) / i
          end do
          pp = (p1 - z * p0) * (nx / (1.0_RD_ - z * z))
          dz = -p0 / pp
          z = z + dz
          if (abs(dz) < RDEPS_) exit
       end do
       x(ih) = -z
       w(ih) = 1.0_RD_ / ((1.0_RD_ - z * z) * pp * pp)
       x(nx - ih + 1) = z
       w(nx - ih + 1) = w(ih)
    end do

  end subroutine gaussLegen


  !+
  ! Interpolate & get a finer gridded vector
  !-
  function grid_refine_1R(dat0, nsub, mfunc) result(dat1)

    real(R_), intent(in) :: dat0(:) ! original tabulated data
    integer,  intent(in) :: nsub    ! # of subintervals
    integer,  intent(in) :: mfunc   ! function flag [0,2]
    real(R_) :: dat1(nsub * size(dat0) - nsub + 1) ! interpolated data
    ! -   mfunc = 0 : linear     ( Y = a * X + b )
    !     mfunc = 1 : log-linear ( Y = a * logX + b, X > 0)
    !     mfunc = 2 : log-log    ( Y = b * exp(a * X), logY = log(b) + a * X, X > 0)
    integer  :: nx0, nx1, isub, ix0, ix1
    real(R_) :: x, x1, x2

    nx0 = size(dat0)
    nx1 = nsub * (nx0 - 1) + 1
    x1 = 0.0_R_
    x2 = real(nsub, R_)
    do ix0 = 1, nx0 - 1
       do isub = 1, nsub
          x = real(isub - 1, R_)
          ix1 = nsub * (ix0 - 1) + isub
          dat1(ix1) = intp_R(mfunc, x1, x2, dat0(ix0), dat0(ix0 + 1), x)
       end do
    end do
    dat1(nx1) = dat0(nx0)

  end function grid_refine_1R


  !+
  ! Integrate y(x)*dx, where dx = h, with trapezoidal rule
  !-
  function integ_trapz_eqX_R(y, h) result(res)

    real(R_), intent(in) :: y(:) ! y data
    real(R_), intent(in) :: h    ! x spacing (constant for all intervals)
    real(R_)  :: res ! integral
    integer :: n

    n = size(y)
    res = h * ((y(1) + y(n)) * 0.5_R_ + sum(y(2:n-1)))

  end function integ_trapz_eqX_R


  !+
  ! Integrate y(x)*dx with trapezoidal rule
  !-
  function integ_trapz_R(x, y) result(res)

    real(R_), intent(in) :: x(:), y(:) ! x & y data
    real(R_)  :: res ! integral
    integer :: n
    n = size(y)
    res = y(1) * (x(2) - x(1)) + y(n) * (x(n) - x(n-1))
    res = 0.5_R_ * (res + sum(y(2:n-1) * (x(3:n) - x(1:n-2))))

  end function integ_trapz_R


  !+
  ! Interpolation with following functions
  !-
  function intp_R(mfunc, x1, x2, y1, y2, x) result(y)

    integer,  intent(in) :: mfunc  ! flag for function type (see below)
    real(R_), intent(in) :: x1, x2 ! two X grid values
    real(R_), intent(in) :: y1, y2 ! two Y grid values
    real(R_), intent(in) :: x      ! X data
    ! -   mfunc = 0 : linear     ( Y = a * X + b )
    !     mfunc = 1 : log-linear ( Y = a * logX + b, X > 0)
    !     mfunc = 2 : log-log    ( Y = b * exp(a * X), logY = log(b) + a * X, X > 0 )
    real(R_) :: y, a, dx

    if (mfunc == 0) then ! linear
       dx = x2 - x1
       if (abs(dx) < RTINY_) then
          y = y1
       else
          a = (x - x1) / dx
          y = (1.0_R_ - a) * y1 + a * y2
       end if

    else
       dx = log(x2) - log(x1)
       if (abs(dx) < RTINY_) then
          y = y1
       else
          a = (log(x) - log(x1)) / dx
          if (mfunc == 1) then ! log-linear
             y = (1.0_R_ - a) * y1 + a * y2
          else ! log-log
             y = (1.0_R_ - a) * log(y1) + a * log(y2)
             y = exp(y)
          end if
       end if
    end if

  end function intp_R


  !+
  ! Akima interpoaltion: To get coefficients (derivatives, dy/dx)
  !
  ! Reference: Akima, H., 1970: J. ACM, 17, 589-602.
  !-
  function intp_Akima_coef_1R(x, y) result(res)

    real(R_), intent(in) :: x(:) ! tabulated x
    real(R_), intent(in) :: y(:) ! tabulated y
    real(R_) :: res(size(y)) ! dy/dx, derivative estimates
    integer  :: n, i
    real(R_) :: d(-1 : size(y)+1), w0, w1

    ! Trapezoidal differences
    n = size(y)
    d(1:n-1) = (y(2:n) - y(1:n-1)) / (x(2:n) - x(1:n-1))
    d(0)   = 2.0_R_ * d(1)   - d(2)
    d(-1)  = 2.0_R_ * d(0)   - d(1)
    d(n)   = 2.0_R_ * d(n-1) - d(n-2)
    d(n+1) = 2.0_R_ * d(n)   - d(n-1)

    ! Derivative estimates
    do i = 1, n
       w1 = abs(d(i+1) - d(i)  )
       w0 = abs(d(i-1) - d(i-2))
       if (w1 + w0 == 0.0_R_) then
          res(i) = (d(i-1) + d(i)) * 0.5_R_
       else
          res(i) = (w1 * d(i-1) + w0 * d(i)) / (w1 + w0)
       end if
    end do

  end function intp_Akima_coef_1R


  !+
  ! Bilinear interpolation for a 2D grid
  !-
  function intp_bilinear_R(t, u, dat) result(res)

    real(R_), intent(in) :: t, u  ! normalized (x, y) coordinates in [0,1]
    !// t = [0,1] for data = [dat(1,:), dat(2,:)]
    !   u = [0,1] for data = [dat(:,1), dat(:,2)]
    !   Note: t & u would be in [0,1] usually. It is implicitly assumed in most cases.
    !   However, there is no such a limitation in the algorithm used in this code.
    real(R_), intent(in) :: dat(:,:) ! data values at 2*2 points
    real(R_) :: res ! result, interpolated data value at (t, u)
    real(R_) :: daty(2)

    daty(1:2) = (1.0_R_ - t) * dat(1, 1:2) + t * dat(2, 1:2)
    res = (1.0_R_ - u) * daty(1) + u * daty(2)

  end function intp_bilinear_R


  !+
  ! Bicubic interpolation for equally-spaced 2D grid
  !-
  function intp_bicubic_eqXY_R(t, u, dat) result(res)

    real(R_), intent(in) :: t, u  ! normalized (x, y) coordinates
    !// t = [-1, 2] for data = [dat(1,:), dat(4,:)]
    !   u = [-1, 2] for data = [dat(:,1), dat(:,4)]
    !   Note: t & u would be in [0,1] usually. It is implicitly assumed in most cases.
    !   However, there is no such a limitation in the algorithm used in this code.
    real(R_), intent(in) :: dat(:,:) ! data values at 4*4 points
    real(R_) :: res ! result, interpolated data value at (t, u)
    real(R_) :: a23, a14, b23, b14, a(4), b(4), daty(4)
    integer :: iy

    ! Cefficients
    a14 = FRAC16_ * (t - 1.0_R_) * t
    b14 = FRAC16_ * (u - 1.0_R_) * u
    a23 = 0.5_R_  * (t - 2.0_R_) * (t + 1.0_R_)
    b23 = 0.5_R_  * (u - 2.0_R_) * (u + 1.0_R_)
    a(1) = -a14 * (t - 2.0_R_)
    a(2) =  a23 * (t - 1.0_R_)
    a(3) = -a23 * t
    a(4) =  a14 * (t + 1.0_R_)
    b(1) = -b14 * (u - 2.0_R_)
    b(2) =  b23 * (u - 1.0_R_)
    b(3) = -b23 * u
    b(4) =  b14 * (u + 1.0_R_)

    ! Interpolate
    do iy = 1, 4
       daty(iy) = dot_product(dat(1:4, iy), a(1:4))
    end do
    res = dot_product(daty, b)

  end function intp_bicubic_eqXY_R


  !+
  ! Interpolation by cubic polynomial, using 3rd-order Lagrangian formula
  ! Note: x-values for the four points should be different each other
  !-
  function intp_cubic_R(x, x4, y4) result(y)

    real(R_), intent(in) :: x             ! x data for the target
    real(R_), intent(in) :: x4(4), y4(4)  ! data for four points
    real(R_) :: y                         ! interpolated y value
    real(R_) :: yy1, yy2, yy3, yy4

    yy1 = y4(1) * (x - x4(2)) / ((x4(1) - x4(2)) * (x4(1) - x4(3)) * (x4(1) - x4(4)))
    yy2 = y4(2) * (x - x4(1)) / ((x4(2) - x4(3)) * (x4(2) - x4(4)) * (x4(2) - x4(1)))
    yy3 = y4(3) * (x - x4(4)) / ((x4(3) - x4(4)) * (x4(3) - x4(1)) * (x4(3) - x4(2)))
    yy4 = y4(4) * (x - x4(3)) / ((x4(4) - x4(1)) * (x4(4) - x4(2)) * (x4(4) - x4(3)))
    y = (x - x4(3)) * (x - x4(4)) * (yy1 + yy2) + (x - x4(1)) * (x - x4(2)) * (yy3 + yy4)

    !real(R_) :: y12, y23, y34, y123, y234
    !y12  = ((x - x4(2)) * y4(1) - (x - x4(1)) * y4(2)) / (x4(1) - x4(2))
    !y23  = ((x - x4(3)) * y4(2) - (x - x4(2)) * y4(3)) / (x4(2) - x4(3))
    !y34  = ((x - x4(4)) * y4(3) - (x - x4(3)) * y4(4)) / (x4(3) - x4(4))
    !y123 = ((x - x4(3)) * y12   - (x - x4(1)) * y23  ) / (x4(1) - x4(3))
    !y234 = ((x - x4(4)) * y23   - (x - x4(2)) * y34  ) / (x4(2) - x4(4))
    !y    = ((x - x4(4)) * y123  - (x - x4(1)) * y234 ) / (x4(1) - x4(4))

  end function intp_cubic_R


  !+
  ! Interpolation by cubic polynomial for equally-spaced x-grids
  ! X should be relatively scaled in [-1,2] for the range of [y(1),y(4)]
  ! Thus, x = [0,1] for y=[y(2),y(3)].
  !-
  function intp_cubic_eqX_R(xrat, y4) result(y)

    real(R_), intent(in) :: xrat  ! relative fraction of x ([-1,2] for y=[y4(1),y4(4)])
    real(R_), intent(in) :: y4(:) ! y values at 4 points
    real(R_) :: y ! result, interpolated y value
    real(R_) :: a23, a14

    a23 = 0.5_R_  * (xrat - 2.0_R_) * (xrat + 1.0_R_)
    a14 = FRAC16_ * (xrat - 1.0_R_) * xrat
    y = a14 * (-(xrat - 2.0_R_) * y4(1) + (xrat + 1.0_R_) * y4(4)) &
         & + a23 * ((xrat - 1.0_R_) * y4(2) - xrat * y4(3))

    !real(R_) :: y12, y23, y34, y123, y234
    !y12  =          - xrat  * y4(1) + (xrat + 1.0_R_) * y4(2) !BUGFIX, 16/01/2009
    !y23  =  (1.0_R_ - xrat) * y4(2) +  xrat           * y4(3)
    !y34  =  (2.0_R_ - xrat) * y4(3) + (xrat - 1.0_R_) * y4(4)
    !y123 =  (1.0_R_ - xrat) * y12   + (xrat + 1.0_R_) * y23
    !y234 =  (2.0_R_ - xrat) * y23   +  xrat           * y34
    !y    = ((2.0_R_ - xrat) * y123  + (xrat + 1.0_R_) * y234) / 6.0_R_

  end function intp_cubic_eqX_R


  !+
  ! Interpolation for a point by 3rd-order Hermite polynomial
  !-
  function intp_Hermite3_R(x0, x1, y0, y1, yp0, yp1, x) result(res)

    real(R_), intent(in) :: x0, x1   ! x values at the bounds
    real(R_), intent(in) :: y0, y1   ! y values at x0 and x1
    real(R_), intent(in) :: yp0, yp1 ! dy/dx values at x0 and x1
    real(R_), intent(in) :: x        ! x value at the target
    !// x would be usually in [x0, x1]. However, there is no explicit limitation in this code.
    real(R_) :: res ! y value interpolated for the x
    real(R_) :: b, c, d, dx, t

    dx = x1 - x0
    b = yp0 * dx
    call intp_Hermite3_coef(y0, y1, b, yp1*dx, c, d)
    t = (x - x0) / dx
    res = y0 + t * (b + t * (c + d * t))

  end function intp_Hermite3_R


  !+
  ! Interpolation for a point by 3rd-order Hermite polynomial
  !// Note: This can be used for a Bezier curve, by scaling the derivatives:
  !     (yp0, yp1) --> (yp0 * fac, yp1 * fac)
  !   The scaling factor (fac) is usually 3.
  !-
  function intp_Hermite3_1R(x0, x1, y0, y1, yp0, yp1, x) result(res)

    real(R_), intent(in) :: x0, x1   ! x values at the bounds
    real(R_), intent(in) :: y0, y1   ! y values at x0 and x1
    real(R_), intent(in) :: yp0, yp1 ! dy/dx values at x0 and x1
    real(R_), intent(in) :: x(:)     ! x values at the target points
    !// x would be usually in [x0, x1]. However, there is no explicit limitation in this code.
    real(R_) :: res(size(x)) ! y values interpolated for the x
    real(R_) :: t(size(x))
    real(R_) :: b, c, d, dx

    dx = x1 - x0
    b = yp0 * dx
    call intp_Hermite3_coef(y0, y1, b, yp1*dx, c, d) ! for t = [0,1]
    t(:) = (x(:) - x0) / dx
    res(:) = y0 + t(:) * (b + t(:) * (c + t(:) * d))

  end function intp_Hermite3_1R


  !+
  ! Returns 3rd-order Hermite polynomial coefficients
  !-
  subroutine intp_Hermite3_coef(y0, y1, yp0, yp1, c, d)

    real(R_), intent(in) :: y0, y1   ! y     values at t = 0 and 1
    real(R_), intent(in) :: yp0, yp1 ! dy/dt values at t = 0 and 1
    real(R_), intent(out) :: c, d ! coefficients
    !// Hermite_polynomial = y0 + yp0 * t + c * t**2 + d * t**3, for t = [0,1]
    c = 3.0_R_ * (y1 - y0) - (yp0 + yp1) - yp0
    d = 2.0_R_ * (y0 - y1) + (yp0 + yp1)

  end subroutine intp_Hermite3_coef


  !+
  ! Linear interpolation from a table with equally-spaced x values
  !-
  function intp_tabLin_eqX_R(ytab, nx, xs, xe, x) result (y)

    real(R_), intent(in) :: ytab(:) ! y tabulated values
    integer,  intent(in) :: nx      ! # of points in table
    real(R_), intent(in) :: xs, xe  ! x values at start and end points [ix = 1 and nx]
    real(R_), intent(in) :: x       ! x value at target point
    real(R_) :: y ! result, interpolated y value
    real(R_) :: rat
    integer  :: ix

    rat = max(0.0_R_, min(1.0_R_, (x - xs) / (xe - xs))) * (nx - 1) ! [0,nx-1)
    ix = min(nx - 2, int(rat)) ! [0,nx-2]
    rat = rat - ix
    y = (1.0_R_ - rat) * ytab(ix + 1) + rat * ytab(ix + 2)

  end function intp_tabLin_eqX_R


  !+
  ! Cubic interpolation from a table with equally-spaced x values
  !-
  function intp_tabCub_eqX_R(ytab, nx, xs, xe, x) result (y)

    real(R_), intent(in) :: ytab(:) ! y tabulated values
    integer,  intent(in) :: nx      ! # of points in table (should be >= 4)
    real(R_), intent(in) :: xs, xe  ! x values at start and end points [ix = 1 and nx]
    real(R_), intent(in) :: x       ! x value at target point
    real(R_) :: y ! result, interpolated y value
    real(R_) :: rat
    integer  :: ix

    rat = max(0.0_R_, min(1.0_R_, (x - xs) / (xe - xs))) * (nx - 1) ! [0,nx-1)
    ix = min(nx - 2, int(rat)) ! [0,nx-2]
    rat = rat - ix
    ix = ix + 1 ! [1,nx-1]
    if (ix >= 2 .and. ix <= nx - 2) then
       y = intp_cubic_eqX_R(rat, ytab(ix-1:ix+2))
    else if (ix <= 1) then
       rat = rat - 1.0_R_
       y = intp_cubic_eqX_R(rat, ytab(1:4))
    else
       rat = rat + 1.0_R_
       y = intp_cubic_eqX_R(rat, ytab(nx-3:nx))
    end if

  end function intp_tabCub_eqX_R


  !+
  ! Interpolation of a scalar by natural cubic spline
  !-
  function intp_spline_R(ax, dx, y0, y1, ddy0, ddy1) result(y)

    real(R_), intent(in)  :: ax         ! ratio, (x - x0) / dx
    real(R_), intent(in)  :: dx         ! x spacing, x1 - x0
    real(R_), intent(in)  :: y0, y1     ! y values at lower and upper points
    real(R_), intent(in)  :: ddy0, ddy1 ! 2nd derivatives at lower and upper points
    real(R_) :: y  ! result, interpolated value of y
    real(R_) :: a
    a = (1.0_R_ - ax) * ax / 6.0_R_ * dx**2
    y = (1.0_R_ - ax) * y0 + ax * y1 - a * ((2.0_R_ - ax) * ddy0 + (1.0_R_ + ax) * ddy1)

  end function intp_spline_R


  !+
  ! Interpolation of a vector by natural cubic spline
  !-
  function intp_spline_1R(ax, dx, y0, y1, ddy0, ddy1) result(y)

    real(R_), intent(in)  :: ax         ! ratio, (x - x0) / dx
    real(R_), intent(in)  :: dx         ! x spacing, x1 - x0
    real(R_), intent(in)  :: y0(:), y1(:)     ! y values at lower and upper points
    real(R_), intent(in)  :: ddy0(:), ddy1(:) ! 2nd derivatives of y at lower and upper points
    real(R_) :: y(size(y0))  ! result, interpolated value of y
    real(R_) :: a, b, c
    a = ax / 6.0_R_ * dx**2
    b = a * (1.0_R_ - ax) * (2.0_R_ - ax)
    c = a * (1.0_R_ - ax**2)
    y(:) = (1.0_R_ - ax) * y0(:) + ax * y1(:) - (b * ddy0(:) + c * ddy1(:))

  end function intp_spline_1R


  !+
  ! Second derivatives of a natural cubic spline
  !-
  function intp_spline_coef_1R(dx, yv) result(y2)

    real(R_), intent(in)  :: dx(:)  ! x spacings, i=1,n-1
    real(R_), intent(in)  :: yv(:)  ! y values,   i=1,n
    real(R_) :: y2(size(yv, 1))  ! 2nd derivatives, d^2(y)/dx^2, i=1,n
    integer  :: i, n
    real(R_) :: p
    real(R_), allocatable :: uv(:)

    ! Initialize
    n = size(yv, 1)
    allocate (uv(n))

    ! Decomposite the tridiagonal matrix
    y2(1) = 0.0_R_
    uv(1) = 0.0_R_
    do i = 2, n - 1
       p = y2(i-1) * dx(i-1) + 2.0_R_ * (dx(i-1) + dx(i))
       y2(i) = -dx(i) / p
       uv(i) = (yv(i+1) - yv(i)) / dx(i) - (yv(i) - yv(i-1)) / dx(i-1)
       uv(i) = (6.0_R_ * uv(i) - dx(i-1) * uv(i-1)) / p
    end do

    ! Backward substitution
    y2(n) = 0.0_R_
    do i = n-1, 1, -1
       y2(i) = y2(i) * y2(i+1) + uv(i)
    end do
    deallocate (uv)

  end function intp_spline_coef_1R


  !+
  ! Interpolation for a point by Lagrange polynominal
  !  This is a numerically intensive method. The # of data points used should not be very large.
  !-
  function intpTab_Lagr_R(x, y, xobj) result(res)

    real(R_), intent(in) :: x(:), y(:) ! x & y data
    real(R_), intent(in) :: xobj ! x at the target point
    real(R_)  :: res
    real(RD_) :: tmp
    real(R_)  :: fl, fli
    integer   :: i, j, n

    n = size(x)
    tmp = 0.0_RD_
    do i = 1, n
       fli = 1.0_R_
       fl  = 1.0_R_
       do j = 1, i - 1
          fli = fli * (xobj - x(j))
          fl  = fl  * (x(i) - x(j))
       end do
       do j = i + 1, n
          fli = fli * (xobj - x(j))
          fl  = fl  * (x(i) - x(j))
       end do
       tmp = tmp + y(i) * (fli / fl)
    end do
    res = tmp

  end function intpTab_Lagr_R


  !+
  ! Akima interpolation for a point, based on a table
  !-
  function intpTab_Akima_R(xtab, ytab, x, mext, iix) result(y)

    real(R_), intent(in) :: xtab(:), ytab(:)  ! X and Y table values
    real(R_), intent(in) :: x                 ! X value at a point of interest
    integer,  optional, intent(in) :: mext    ! if 1, extraporation
    integer,  optional, intent(in) :: iix     ! if present, iix is known X grid index (ix),
    !--                                         and the above "x" should be normalized value in the interval [x(ix),x(ix)]
    real(R_) :: y ! Y values interpolated
    call intpTab_Akima_oR(xtab, ytab, x, y, mext=mext, iix=iix)

  end function intpTab_Akima_R


  !+
  ! Akima interpolation for a point, based on a table
  !
  ! Note: The iix is useful when there are many ytab's for a single set of xtab and x:
  !   do i = 1, size(x)
  !      call gridIdx_loc(xtab, x(i), iix(i), ax(i), mext)
  !   end do
  !   call intpTab_Akima_o1R(xtab, ytab1, x, y1, mext, iix, dydx)
  !   call intpTab_Akima_o1R(xtab, ytab2, x, y2, mext, iix, dydx)
  !   call intpTab_Akima_o1R(xtab, ytab3, x, y3, mext, iix, dydx)
  !   ...
  !-
  subroutine intpTab_Akima_oR(xtab, ytab, x, y, mext, iix, dydx)

    real(R_), intent(in) :: xtab(:), ytab(:) ! X and Y table values
    real(R_), intent(in) :: x                ! X value at the point of interest
    real(R_), intent(out):: y                ! Y value interpolated at X
    integer,  optional, intent(in) :: mext   ! if 1, extraporation
    integer,  optional, intent(in) :: iix    ! if present, iix is known X grid index (ix),
    !--                                        and the above "x" should be normalized value in the interval [x(ix),x(ix)]
    real(R_), optional, intent(out):: dydx   ! dY/dX at X's
    real(R_) :: dydxtab(0:1), d(-2:2), ax, dx, c1, c2, c3, w0, w1
    integer  :: ix, i, nx

    ! Locate in a grid
    nx = size(ytab)
    if (present(iix)) then
       ix = iix
       ax = x
    else
       call gridIdx_loc(xtab, x, ix, ax, mext)
    end if

    ! Trapezoidal differences
    if (ix >= 3 .and. ix <= nx - 3) then
       d(-2:2) = (ytab(ix-1:ix+3) - ytab(ix-2:ix+2)) / (xtab(ix-1:ix+3) - xtab(ix-2:ix+2))
    else if (ix <= 2) then
       if (ix == 2) then
          d(-1:2) = (ytab(ix  :ix+3) - ytab(ix-1:ix+2)) / (xtab(ix  :ix+3) - xtab(ix-1:ix+2))
          d(-2) = 2.0_R_ * d(-1) - d(0)
       else
          d(0:2)  = (ytab(ix+1:ix+3) - ytab(ix  :ix+2)) / (xtab(ix+1:ix+3) - xtab(ix  :ix+2))
          d(-1) = 2.0_R_ * d(0)  - d(1)
          d(-2) = 2.0_R_ * d(-1) - d(0)
       end if
    else
       if (ix == nx - 2) then
          d(-2:1) = (ytab(ix-1:ix+2) - ytab(ix-2:ix+1)) / (xtab(ix-1:ix+2) - xtab(ix-2:ix+1))
          d(2) = 2.0_R_ * d(1) - d(0)
       else !if (ix == nx - 1) then
          d(-2:0) = (ytab(ix-1:ix+1) - ytab(ix-2:ix  )) / (xtab(ix-1:ix+1) - xtab(ix-2:ix  ))
          d(1) = 2.0_R_ * d(0) - d(-1)  
          d(2) = 2.0_R_ * d(1) - d(0)
       end if
    end if

    ! Derivative estimates
    do i = 0, 1
       w1 = abs(d(i+1) - d(i  ))
       w0 = abs(d(i-1) - d(i-2))
       if (w1 + w0 == 0.0_R_) then
          dydxtab(i) = (d(i-1) + d(i)) * 0.5_R_
       else
          dydxtab(i) = (w1 * d(i-1) + w0 * d(i)) / (w1 + w0)
       end if
    end do

    ! Cubic polynomial
    dx = xtab(ix+1) - xtab(ix)
    c1 = dydxtab(0) * dx
    call intp_Hermite3_coef(ytab(ix), ytab(ix+1), c1, dydxtab(1)*dx, c2, c3) ! get coeffs.
    y = ytab(ix) + ax * (c1 + ax * (c2 + ax * c3))
    if (present(dydx)) dydx = (c1 + ax * (2.0_R_ * c2 + 3.0_R_ * c3 * ax)) / dx

  end subroutine intpTab_Akima_oR



  !+
  ! Akima interpolation for points, based on a table
  !-
  function intpTab_Akima_1R(xtab, ytab, x, mext, iix) result(y)

    real(R_), intent(in) :: xtab(:), ytab(:)  ! X and Y table values
    real(R_), intent(in) :: x(:)              ! X values at the points of interest
    integer,  optional, intent(in) :: mext    ! if 1, extraporation
    integer,  optional, intent(in) :: iix(:)  ! if present, iix is known X grid index (ix),
    !--                                         and the above "x" should be normalized value in the interval [x(ix),x(ix)]
    real(R_) :: y(size(x))                    ! Y values interpolated
    call intpTab_Akima_o1R(xtab, ytab, x, y, mext, iix)

  end function intpTab_Akima_1R


  !+
  ! Akima interpolation for points, based on a table
  !
  ! Note: The iix is useful when there are many ytab's for a single set of xtab and x:
  !   do i = 1, size(x)
  !      call gridIdx_loc(xtab, x(i), iix(i), ax(i), mext)
  !   end do
  !   call intpTab_Akima_o1R(xtab, ytab1, x, y1, mext, iix, dydx)
  !   call intpTab_Akima_o1R(xtab, ytab2, x, y2, mext, iix, dydx)
  !   call intpTab_Akima_o1R(xtab, ytab3, x, y3, mext, iix, dydx)
  !   ...
  !-
  subroutine intpTab_Akima_o1R(xtab, ytab, x, y, mext, iix, dydx)

    real(R_), intent(in) :: xtab(:), ytab(:)  ! X and Y table values
    real(R_), intent(in) :: x(:)              ! X values at the points of interest
    real(R_), intent(out):: y(:)              ! Y values interpolated at X's
    integer,  optional, intent(in) :: mext    ! if 1, extraporation
    integer,  optional, intent(in) :: iix(:)  ! if present, iix is known X grid index (ix),
    !--                                         and the above "x" should be normalized value in the interval [x(ix),x(ix)]
    real(R_), optional, intent(out):: dydx(:) ! dY/dX at X's
    real(R_) :: dydxtab(size(xtab)), ax, dx, b, c, d
    integer  :: ix, i

    dydxtab(:) = intp_Akima_coef_1R(xtab, ytab) ! derivatives    
    do i = 1, size(x)

       ! Grid index and interpolation factor
       if (present(iix)) then
          ix = iix(i)
          ax = x(i)
       else
          call gridIdx_loc(xtab, x(i), ix, ax, mext) ! locate in a grid
       end if

       ! Interpolate Y(X) and possibly get dY/dX at X
       dx = xtab(ix+1) - xtab(ix)
       b = dydxtab(ix) * dx
       call intp_Hermite3_coef(ytab(ix), ytab(ix+1), b, dydxtab(ix+1)*dx, c, d) ! cubic polynomial coeffs.
       y(i) = ytab(ix) + ax * (b + ax * (c + ax * d))
       if (present(dydx)) dydx(i) = (b + ax * (2.0_R_ * c + 3.0_R_ * d * ax)) / dx
    end do

  end subroutine intpTab_Akima_o1R


  !+
  ! Linear interpolation for a point, based on a table
  !-
  function intpTab_linear_R(xtab, ytab, x, mext) result(y)

    real(R_), intent(in) :: xtab(:), ytab(:) ! X and Y table values
    real(R_), intent(in) :: x ! X value at the point of interest
    integer,  intent(in), optional :: mext ! if 1, extraporation
    real(R_) :: y ! Y value interpolated
    real(R_) :: rat
    integer  :: i
    call gridIdx_loc(xtab, x, i, rat, mext)
    y = (1.0_R_ - rat) * ytab(i)  + rat * ytab(i + 1)

  end function intpTab_linear_R


  !+
  ! Linear interpolation for points, based on a table
  !-
  function intpTab_linear_1R(xtab, ytab, x, mext) result(y)

    real(R_), intent(in) :: xtab(:), ytab(:) ! X and Y table values
    real(R_), intent(in) :: x(:) ! X values at the points of interest
    integer,  intent(in), optional :: mext ! if 1, extraporation
    real(R_) :: y(size(x)) ! Y values interpolated
    real(R_) :: rat
    integer  :: i, idat
    do idat = 1, size(x)
       call gridIdx_loc(xtab, x(idat), i, rat, mext)
       y(idat) = (1.0_R_ - rat) * ytab(i)  + rat * ytab(i + 1)
    end do

  end function intpTab_linear_1R


  !+
  ! Cubic spline interpolation for points, based on a table
  !-
  function intpTab_spline_1R(xtab, ytab, x, mext) result(y)

    real(R_), intent(in) :: xtab(:), ytab(:) ! X and Y table values
    real(R_), intent(in) :: x(:) ! X values at the points of interest
    integer,  intent(in), optional :: mext ! if 1, extraporation
    real(R_) :: y(size(x)) ! Y values interpolated
    real(R_) :: ax, dx(size(xtab)-1), dy2(size(xtab)), a
    integer  :: ix, i, nx

    ! Initialize
    nx = size(ytab)
    dx(1:nx-1) = xtab(2:nx) - xtab(1:nx-1)
    dy2(:) = intp_spline_coef_1R(dx, ytab)

    ! Interpolate for each X point
    do i = 1, size(x)
       call gridIdx_loc(xtab, x(i), ix, ax, mext)
       a = (1.0_R_ - ax) * ax / 6.0_R_ * dx(ix)**2
       y(i) = (1.0_R_ - ax) * ytab(ix) + ax * ytab(ix+1) &
            - a * ((2.0_R_ - ax) * dy2(ix) + (1.0_R_ + ax) * dy2(ix+1))
    end do

  end function intpTab_spline_1R


  !+
  ! Legendre polynomials for i=[1:n]
  !-
  function legenPoly_1R(x, n) result(pn)

    real(R_), intent(in) :: x  ! X in [-1, 1]
    integer,  intent(in) :: n  ! max order (should be >= 1)
    real(R_) :: pn(n)  ! Legendre polynomials for orders i=[1:n]
    ! -  Note: P_n = 1 for n=0
    real(R_) :: p0, p1
    integer  :: i, ii

    ! The first terms
    p0 = 1.0_R_
    p1 = x
    pn(1) = x

    ! Higher terms
    ii = 1
    do i = 2, n
       ii = ii + 2
       pn(i) = (real(ii, R_) * x * p1 - real(i-1, R_) * p0) / real(i, R_)
       p0 = p1
       p1 = pn(i)
    end do

  end function legenPoly_1R


  !+
  ! Form an expansion with Legendre polynomials from expansion coefficients
  !-
  function legenPolySer_1R(g0, gn, x) result(p)

    real(R_), intent(in) :: g0       ! 0th moment
    real(R_), intent(in) :: gn(:)    ! moments (1st and higher)
    real(R_), intent(in) :: x(:)     ! tabulated x values in [-1,1]
    real(R_) :: p(size(x))  ! result
    real(R_) :: p0(size(x)), p1(size(x)), p2(size(x))
    integer :: ig, ng, i2

    p(:) = g0
    ng = size(gn)
    if (ng >= 1) then
       p0(:) = 1.0_R_
       p1(:) = x(:)
       i2 = 3
       p(:) = p(:) + real(i2, R_) * gn(1) * p1(:)
       do ig = 2, ng
          p2(:) = (real(i2, R_) * x(:) * p1 - real(ig-1, R_) * p0(:)) / real(ig, R_)
          i2 = i2 + 2
          p(:) = p(:) + real(i2, R_) * gn(ig) * p2(:)
          p0(:) = p1(:)
          p1(:) = p2(:)
       end do
    end if

  end function legenPolySer_1R


  !+
  ! Modify integral of y*dx as prescribed
  !-
  function modifyTotal_1R(x, y, tot) result(y1) 

    real(R_), intent(in) :: x(:), y(:) ! tabulated (x,y) data
    real(R_), intent(in) :: tot ! prescribed integral of y(x)*dx 
    real(R_) :: y1(size(y))
    y1(:) = y(:) * (integ_trapz_R(x, y) / tot)

  end function modifyTotal_1R


  !+
  ! Vector of relative deviation (from 1) increased by a small forcing
  !    y = x + (1 + x) * f
  !  where y << 1, x << 1, f << 1.
  ! Purpose: to accurately compute glowth of relative deviation due to (very) small forcing
  !       Y =      X  * (1 + f), where Y = 1 + y, X = 1 + x, so that
  !   1 + y = (1 + x) * (1 + f).
  !  Because the above is not safe to compute accurate y, we use instead
  !       y = x + (1 + x) * f
  !-
  function relaDevMod_1R(x, f) result(y)

    real(R_), intent(in) :: x(:) ! original relative deviation from 1
    real(R_), intent(in) :: f(:) ! forcing
    real(R_) :: y(size(x, 1))    ! modified relative deviation from 1
    y(:) = x(:) + (1.0_R_ + x(:)) * f(:)

  end function relaDevMod_1R


  !+
  ! Solve a quadratic equation
  !  Y = a * X^2 + b * X + c = 0
  !-
  subroutine root_quad(a, b, c, s1, s2, dd)

    real(R_), intent(in)  :: a, b, c ! A, B, C
    real(R_), intent(out) :: s1  ! one solution
    real(R_), intent(out) :: s2  ! the other solution, abs(s2) <= abs(s1)
    real(R_), optional, intent(out) :: dd  ! determinant D
    real(R_) :: d

    d = b**2 - 4.0_R_ * a * c
    if (abs(a) <= abs(b) * RSML_) then ! A ~= 0
       s1 = -c / b
       s2 = s1
    else
       if (d > 0.0_R_) then
          s1 = 0.5_R_ * (-b + sign(sqrt(d), -b)) / a
          s2 = c / (a * s1)
       else if (d == 0.0_R_) then
          s1 = -0.5_R_ * b / a
          s2 = s1
       else
          s1 = 0.0_R_ ! zero, tentatively
          s2 = 0.0_R_
       end if
    end if
    if (present(dd)) dd = d

  end subroutine root_quad


  !+
  ! Solve a quadratic equation with the first coefficient of unity (A=1)
  !  Y = X^2 + b * X + c = 0
  !-
  subroutine root_quadA1(b, c, s1, s2)

    real(R_), intent(in)  :: b   ! B
    real(R_), intent(in)  :: c   ! C
    real(R_), intent(out) :: s1  ! one solution
    real(R_), intent(out) :: s2  ! the other solution, abs(s2) <= abs(s1)
    real(R_) :: d

    d = b**2 - 4.0_R_ * c
    if (d > 0.0_R_) then
       s1 = 0.5_R_ * (-b + sign(sqrt(d), -b))
       s2 = c / s1
    else
       s1 = -0.5_R_ * b ! solutions for D=0, tentatively
       s2 = s1 ! if this is the case, it means that there is no solution!
    end if

  end subroutine root_quadA1


  !+
  ! Compute a single root of quadratic equation
  !  Y = X^2 + B * X + C
  !  Root: X = -B (+-) sqrt(B**2 + C)
  !-
  function root_quadA1_one_R(b, c, f) result(res)

    real(R_) :: res             ! -B (+-) sqrt(B**2 + C)
    real(R_), intent(in) :: b   ! B
    real(R_), intent(in) :: c   ! C
    real(R_), intent(in) :: f   ! the sign (+1 or -1)
    real(R_) :: d

    d = b**2 + c
    if (d <= 0.0_R_) then ! invalid input!
       res = -b ! which means there is no solution
    else if (f * b > 0.0_R_) then ! difficult case
       res = -b - sign(sqrt(d), f) 
       res = -c / res
    else                 ! normal case
       res = -b + sign(sqrt(d), f)
    end if

  end function root_quadA1_one_R


  !+
  ! Sectoinal averages for exponential sub-section profile data
  !-
  function secAve_exp_1R(x, po, xa) result(res)

    real(R_), intent(in) :: x(:)  ! original X data
    real(R_), intent(in) :: po(:) ! original profile data (should be > 0)
    real(R_), intent(in) :: xa(:) ! X values at section boundaries for averaged profile
    real(R_) :: res(size(xa) - 1)
    real(R_) :: pp(size(po))
    real(R_) :: rat, ppb, ppt, taua
    integer :: iz, nz, ixb, ixt, ix

    ! Initialize
    nz = size(xa) - 1
    pp(:) = log(max(RTINY_, po(:)))
    call gridIdx_loc(x, xa(1), ixb, rat)
    ppb = (1.0_R_ - rat) * pp(ixb)  + rat * pp(ixb + 1)

    ! Loop for every layers
    do iz = 1, nz
       call gridIdx_loc(x, xa(iz + 1), ixt, rat)
       ppt = (1.0_R_ - rat) * pp(ixt)  + rat * pp(ixt + 1)
       if (ixb == ixt) then
          res(iz) = effAve_exp_R(exp(ppb), exp(ppt))
       else
          taua = 0.0_R_
          taua = taua + effAve_exp_R(exp(ppb), po(ixb + 1)) * (x(ixb + 1) - xa(iz)) ! bottom
          taua = taua + effAve_exp_R(exp(ppt), po(ixt))     * (xa(iz + 1) - x(ixt)) ! top
          do ix = ixb + 1, ixt - 1
             taua = taua + effAve_exp_R(po(ix), po(ix + 1)) * (x(ix + 1) - x(ix))
          end do
          res(iz) = taua / (xa(iz + 1) - xa(iz))
       end if
       ixb  = ixt
       ppb = ppt
    end do

  end function secAve_exp_1R


  !+
  ! Sectoinal averages for linear sub-section profile data
  !-
  function secAve_lin_1R(x, po, xa) result(res)

    real(R_), intent(in) :: x(:)  ! original X data
    real(R_), intent(in) :: po(:) ! original profile data (should be > 0)
    real(R_), intent(in) :: xa(:) ! X values at section boundaries for averaged profile
    real(R_) :: res(size(xa) - 1)
    real(R_) :: rat, ppb, ppt, taua
    integer :: iz, nz, ixb, ixt, ix

    ! Initialize
    nz = size(xa) - 1
    call gridIdx_loc(x, xa(1), ixb, rat)
    ppb = (1.0_R_ - rat) * po(ixb)  + rat * po(ixb + 1)

    ! Loop for every layers
    do iz = 1, nz
       call gridIdx_loc(x, xa(iz + 1), ixt, rat)
       ppt = (1.0_R_ - rat) * po(ixt)  + rat * po(ixt + 1)
       if (ixb == ixt) then
          res(iz) = 0.5_R_ * (ppb + ppt)
       else
          taua = 0.0_R_
          taua = taua + 0.5_R_ * (ppb + po(ixb + 1)) * (x(ixb + 1) - xa(iz)) ! bottom
          taua = taua + 0.5_R_ * (ppt + po(ixt))     * (xa(iz + 1) - x(ixt)) ! top
          do ix = ixb + 1, ixt - 1
             taua = taua + 0.5_R_ * (po(ix) + po(ix + 1)) * (x(ix + 1) - x(ix))
          end do
          res(iz) = taua / (xa(iz + 1) - xa(iz))
       end if
       ixb  = ixt
       ppb = ppt
    end do

  end function secAve_lin_1R


  !+
  ! Smoothing a vector by trianglar filter
  !-
  function smooth_triang_1R(datin, nw) result(datout)

    real(R_), intent(in)  :: datin(:) ! original vector
    integer,  intent(in)  :: nw       ! filter half width 
    !// Full width will be (2*nw - 1). Effective width is just the nw.
    real(R_) :: datout(size(datin)) ! smoothed vector
    integer  :: n, i, i1, i2, iw, j1, j2
    real(RD_) :: sum

    n = size(datin)
    do i = 1, n
       i1 = i
       i2 = i
       sum = datin(i) * nw
       do iw = 1, nw - 1
          i1 = i1 - 1
          i2 = i2 + 1
          j1 = i1
          j2 = i2
          if (i1 <= 0    ) j1 = 2 - i1
          if (i2 >= n + 1) j2 = n - (i2 - n)
          sum = sum + (datin(j1) + datin(j2)) * (nw - iw)
       end do
       datout(i) = sum / real(nw * nw, R_)
    end do

  end function smooth_triang_1R


  !+
  ! Modify a vector by smoothing hot spots
  !  that is, by truncating the original hot spots (large values) and redistribute 
  !  the sum of the truncated parts to all elements
  !-
  subroutine smooth_noHotSpots(mhot, fhot, vec, vhot)

    integer,  intent(in)    :: mhot   ! method flag (0/1) for determining hot spot value
    real(R_), intent(in)    :: fhot   ! factor (mhot=1) or hot spot threshold (mhot=0)
    real(R_), intent(inout) :: vec(:) ! I/O vector
    real(R_), intent(out), optional :: vhot ! hot spot threshold
    !- If mhot = 1, fhot is a factor to average of active (non-zero) data
    !- Else, fhot is hot spot threshold itself.
    !- Limitation: All elements of vec should have the same sign (+/-).
    real(R_) :: ssum, vhot1, vfac

    ! Sum
    ssum = sum(vec)
    if (ssum == 0.0_R_) return ! unexpected condition!

    ! Determine a hot spot threshold
    if (mhot == 0) then
       vhot1 = fhot
    else if (ssum > 0.0_R_) then
       vhot1 = fhot * ssum / count(vec > 0.0_R_)
    else
       vhot1 = fhot * ssum / count(vec < 0.0_R_)
    end if
    if (present(vhot)) vhot = vhot1

    ! Truncate hot spots
    if (vhot1 > 0.0_R_) then
       if (maxval(vec) <= vhot1) return
       vec(:) = min(vhot1, vec(:))
    else
       if (minval(vec) >= vhot1) return
       vec(:) = max(vhot1, vec(:))
    end if

    ! Redistribute
    vfac = ssum / sum(vec)
    vec(:) = vec(:) * vfac

  end subroutine smooth_noHotSpots


  !+
  ! Spherical bessel functions of first kind for a real argument
  !   using the downward recurrence relation
  !-
  function sphBesJ_1RD(nmax, x, mthd) result(res)

    integer,   intent(in) :: nmax ! max order
    real(RD_), intent(in) :: x    ! x
    integer,   intent(in), optional :: mthd ! computation method (0 or 1)
    real(RD_) :: res(0:nmax) ! J_n(x), n = 0,1,2,...,nmax
    integer :: i, nm, m
    real(RD_) :: t1 = 0.0_RD_, t2, t3, f
    real(RD_), parameter :: big = 1.0e+250_RD_, abig = 1.0e-250_RD_

    ! Invalid argument X
    if (x <= 0.0_RD_) then
       res(:) = 0.0_RD_

       ! Very small x
    else if (x <= 1.0e-6_RD_) then
       res(0) = 1.0_RD_
       do i = 1, nmax
          res(i) = res(i - 1) * x / (2 * i + 1)
       end do

       ! Large x
    else
       m = 0
       if (present(mthd)) m = mthd
       if (m == 0) then ! Large x, using classical algorithm
          nm = max(nmax, int(x)) + max(6, 2 + int(4.0_RD_ * x**0.3333333_RD_)) ! starting order
          t3 = 0.0_RD_
          t2 = tiny(1.0_RD_) * 10.0_RD_
          do i = nm, 0, -1 ! downward recurrence
             t1 = (2 * i + 3) / x * t2 - t3
             if (i <= nmax) then
                res(i) = t1
                if (abs(t1) >= big) res(i:nmax) = res(i:nmax) * abig
             end if
             if (abs(t1) >= big) then
                t1 = t1 * abig
                t2 = t2 * abig
             end if
             t3 = t2
             t2 = t1
          end do
          res(0) = sin(x) / x
          t1 = res(0) / t1
          res(1:nmax) = res(1:nmax) * t1 ! rescale

          ! Larger x, using logarithmic derivatives
       else
          nm = max(nmax, int(1.2_RD_ * x)) + 20
          f = 0.0_RD_
          do i = nm - 1, 1, -1 ! downward recurrence
             f = (i + 1) / x - 1.0_RD_ / (f + (i + 1) / x) ! logarithmic derivative of j_i(x)
             if (i <= nmax) res(i) = f ! res() is used as temporary storage
          end do
          res(0) = sin(x) / x ! j_0
          do i = 1, nmax
             res(i) = res(i - 1) / (res(i) + i / x) ! j_i
          end do
       end if
    end if

  end function sphBesJ_1RD


  !+
  ! Spherical bessel functions of second kind for a real argument
  !   using the upward recurrence relation
  !-
  function sphBesY_1RD(nmax, x) result(res)

    integer,   intent(in) :: nmax ! max order
    real(RD_), intent(in) :: x    ! x (should be > 0)
    real(RD_) :: res(0:nmax) ! Y_n(x), n = 0,1,2,...,nmax
    integer :: i
    real(RD_) :: y0, y1, y2, ymin

    ! Invalid argument X
    if (x <= 0.0_RD_) then
       res(:) = -RDHUGE_

       ! Very small x
    else if (x <= 1.0e-6_RD_) then
       ymin = -RDHUGE_ * x
       res(0) = -1.0_RD_ / x
       if (res(0) < ymin) then ! the next will be overflow
          res(1:nmax) = -RDHUGE_
       else
          do i = 1, nmax
             res(i) = res(i - 1) * (2 * i - 1) / x
             if (res(i) < ymin / (2 * i + 1)) then ! the next will be overflow
                res(i+1:nmax) = -RDHUGE_
                exit
             end if
          end do
       end if

       ! Larger x
    else
       y0 = -cos(x) / x
       y1 = (y0 - sin(x)) / x
       res(0) = y0
       res(1) = y1
       do i = 2, nmax ! ascending order
          y2 = (2 * i - 1) * y1 / x - y0
          if (y2 > -RDHUGE_) then
             res(i) = y2
          else
             res(i:nmax) = -RDHUGE_
             exit
          end if
          y0 = y1
          y1 = y2
       end do
    end if

  end function sphBesY_1RD

end module hparx_math
