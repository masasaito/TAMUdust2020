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
! Library of vector and matrix operator procedures
!-
module hparx_vecmat

  use globals
  implicit none
  private

  ! Public
  public :: gridIdx_bin_I     ! find a grid point by binary search
  public :: gridIdx_bin0_I
  public :: gridIdx_bound_I
  public :: gridIdx_fix0
  public :: gridIdx_loc
  public :: gridIdx_unif0_I
  public :: gridIdx_seq_I
  public :: gridIdx_seq0_I
  public :: gridValues_1I
  public :: gridValues_1R
  public :: gridValues_log_1R
  public :: mat_ave_2R
  public :: mat_diag_1R
  public :: mat_diag_2R
  public :: mat_diagProd_R
  public :: mat_diag_scale
  public :: mat_nondiag_scale
  public :: mat_idt_2R
  public :: mat_trace_R
  public :: vec_ave_1R
  public :: vec_killSmall
  public :: vec_norm_R
  public :: vec_outerSum_2R
  public :: vec_outerProd_2R
  public :: vec_outerDiv_2R

contains

  !+
  ! Search the section by the binary search method.
  !  grd(i) <= dat < grd(i+1) or grd(i) >= dat > grd(i+1).
  !   where jmin <= i <= jmax-1
  !  i=jmin   if dat < grd(jmin) or dat >= grd(jmin)
  !  i=jmax-1 if dat > grd(jmax) or dat <= grd(jmax)
  !-
  function gridIdx_bin_I(grd, jmin, jmax, dat) result(i0)

    integer,  intent(in) :: jmin, jmax  ! min & max of the grid index
    real(R_), intent(in) :: grd(:)      ! values at grid points
    real(R_), intent(in) :: dat         ! the target data value 
    integer  :: i0  ! found index, where jmin <= i0 <= jmax-1
    integer  :: i1, i

    i0 = jmin
    i1 = jmax

    if (grd(i0) < grd(i1)) then ! increasing with i
       do
          if (i1 <= i0 + 1) exit
          i = (i0 + i1) / 2
          if (dat >= grd(i)) then
             i0 = i
          else
             i1 = i
          end if
       end do

    else                      ! decreasing with i
       do
          if (i1 <= i0 + 1) exit
          i = (i0 + i1) / 2
          if (dat <= grd(i)) then
             i0 = i
          else
             i1 = i
          end if
       end do
    end if

  end function gridIdx_bin_I


  !+
  ! Binary search a section in grid
  !-
  function gridIdx_bin0_I(grd, dat, iimin, iimax) result(ii)

    real(R_), intent(in) :: grd(0:) ! grid point values, increasing vector
    !// Here, grd(i) should increase with i, and i in [iimin, iimax-1] could be accessed.
    real(R_), intent(in) :: dat     ! target data
    integer,  intent(in) :: iimin   ! min index for ii (>= 0)
    integer,  intent(in) :: iimax   ! max index for ii (<= ubound(grd)+1)
    integer   :: ii  ! result, in [iimin, iimax] 
    !// For ii, grd(ii-1) <= dat < grd(ii), usually. Exceptions are too small/large cases.
    integer   :: i, iu

    if (dat < grd(iimin)) then ! too small
       ii = iimin
    else if (dat >= grd(iimax - 1)) then ! too large
       ii = iimax
    else ! usual case
       ii = iimin + 1
       iu = iimax - 1
       do
          if (iu <= ii) exit
          i = (ii + iu) / 2
          if (dat >= grd(i)) then
             ii = i + 1
          else
             iu = i
          end if
       end do
    end if

  end function gridIdx_bin0_I


  !+
  ! Set an index in a range with open/cyclic/reflective boundary condition
  !-
  function gridIdx_bound_I(mbnd, ix, ixmin, ixmax) result(ix1)

    integer, intent(in) :: mbnd ! flag for boundary condition
    integer, intent(in) :: ix   ! a new index, possibly out of bound
    integer, intent(in) :: ixmin, ixmax ! index min & max
    integer :: ix1

    ix1 = ix

    ! Open
    if (mbnd == 0) then
       if (ix1 > ixmax) ix1 = ixmax + 1
       if (ix1 < ixmin) ix1 = ixmin - 1
       ! Cyclic
    else if (mbnd == 1) then
       if (ix1 > ixmax) ix1 = ix1 - (ixmax - ixmin + 1)
       if (ix1 < ixmin) ix1 = ix1 + (ixmax - ixmin + 1)
       ! Reflective
    else
       if (ix1 > ixmax) ix1 = ixmax - ix1 + 1 + ixmax
       if (ix1 < ixmin) ix1 = ixmin - ix1
    end if

  end function gridIdx_bound_I


  !+
  ! Check an index of grid and possibly fix it
  !-
  subroutine gridIdx_fix0(xgrd, x, ixmax, ix)

    real(R_), intent(in) :: xgrd(0:) ! grid values, increasing vector
    real(R_), intent(in) :: x        ! target data
    integer,  intent(in) :: ixmax    ! max grid index
    integer,  intent(inout) :: ix    ! grid index

    if (x < xgrd(ix - 1)) then
       ix = ix - 1
    else if (x >= xgrd(ix)) then
       ix = min(ix + 1, ixmax)
    end if

  end subroutine gridIdx_fix0


  !+
  ! Locate a point in grid
  !-
  subroutine gridIdx_loc(xg, x, ix, rat, mext)

    real(R_), intent(in)  :: xg(:) ! vector
    real(R_), intent(in)  :: x     ! a value
    integer,  intent(out) :: ix    ! found grid index (1 <= ix <= size(xg)-1)
    real(R_), intent(out) :: rat   ! location between xg(ix) and xg(ix+1)
    !// xg(ix) <= x < xg(ix+1) if xg is  increasing vector and xg(1) <= x < xg(nx)
    !// xg(ix) >= x > xg(ix+1) if xg is deccreasing vector and xg(1) >= x > xg(nx)
    integer,  intent(in), optional :: mext ! if 1, extraporation (rat will be < 0 or > 1)
    integer  :: i, iz

    ! Too few grid points
    iz = size(xg)
    if (iz <= 1) then
       ix = 1
       rat = 0.0_R_

       ! Increasing vector
    else if (xg(1) < xg(iz)) then
       if (x <= xg(1)) then ! out of lower limit
          ix = 1
          rat = 0.0_R_
          if (present(mext)) then
             if (mext == 1) rat = (x - xg(1)) / (xg(2) - xg(1))
          end if
       else if (x >= xg(iz)) then ! out of upper limit
          ix = iz - 1
          rat = 1.0_R_
          if (present(mext)) then
             if (mext == 1) rat = (x - xg(ix)) / (xg(iz) - xg(ix))
          end if
       else                ! in the grid
          ix = 1
          do
             if (iz <= ix + 1) exit
             i = (ix + iz) / 2
             if (x >= xg(i)) then
                ix = i
             else
                iz = i
             end if
          end do
          rat = (x - xg(ix)) / (xg(ix + 1) - xg(ix))
       end if

       ! Decreasing vector
    else
       if (x >= xg(1)) then ! out of lower limit
          ix = 1
          rat = 0.0_R_
          if (present(mext)) then
             if (mext == 1) rat = (x - xg(1)) / (xg(2) - xg(1))
          end if
       else if (x <= xg(iz)) then ! out of upper limit
          ix = iz - 1
          rat = 1.0_R_
          if (present(mext)) then
             if (mext == 1) rat = (x - xg(ix)) / (xg(iz) - xg(ix))
          end if
       else                ! in the grid
          ix = 1
          do
             if (iz <= ix + 1) exit
             i = (ix + iz) / 2
             if (x <= xg(i)) then
                ix = i
             else
                iz = i
             end if
          end do
          rat = (x - xg(ix)) / (xg(ix + 1) - xg(ix))
       end if
    end if

  end subroutine gridIdx_loc


  !+
  ! Sequencial search a section in grid, with an initial estimate
  !    xgrd(ix) <= xdat <  xgrd(ix+1) when xgrd is increasing vector, or
  !    xgrd(ix) >  xdat >= xgrd(ix+1) when xgrd is decreasing vector,
  ! Output ix will be in the range [ixmin, ixmax-1]
  !   ix = ixmin   for xdat <  xgrd(ixmin) or xdat >= xgrd(ixmin) 
  !   ix = ixmax-1 for xdat >= xgrd(ixmax) or xdat <  xgrd(ixmax)
  !-
  function gridIdx_seq_I(xgrd, xdat, ix0, ixmin, ixmax) result(ix)

    real(R_), intent(in) :: xgrd(:) ! grid values
    real(R_), intent(in) :: xdat    ! target data
    integer,  intent(in) :: ix0     ! initial estimate
    integer,  intent(in) :: ixmin   ! min index
    integer,  intent(in) :: ixmax   ! max index
    integer   :: ix ! result, found index

    ! Increasing vector
    if (xgrd(ixmin) < xgrd(ixmax)) then
       if (xdat < xgrd(ix0)) then
          do ix = ix0 - 1, ixmin, -1
             if (xdat >= xgrd(ix)) exit
          end do
          ix = max(ixmin, ix)
       else
          do ix = ix0 + 1, ixmax
             if (xdat < xgrd(ix)) exit
          end do
          ix = min(ixmax - 1, ix - 1)
       end if

       ! Decreasing vector
    else
       if (xdat >= xgrd(ix0)) then
          do ix = ix0 - 1, ixmin, -1
             if (xdat < xgrd(ix)) exit
          end do
          ix = max(ixmin, ix)
       else
          do ix = ix0 + 1, ixmax
             if (xdat >= xgrd(ix)) exit
          end do
          ix = min(ixmax - 1, ix - 1)
       end if
    end if

  end function gridIdx_seq_I


  !+
  ! Sequencial search a section in grid, with an initial estimate
  !-
  function gridIdx_seq0_I(grd, dat, i0, iimin, iimax) result(ii)

    real(R_), intent(in) :: grd(0:) ! grid point values, increasing vector
    !// Here, grd(i) should increase with i, and i in [iimin, iimax-1] could be accessed.
    real(R_), intent(in) :: dat     ! target data
    integer,  intent(in) :: i0      ! initial estimate for ii (1 <= i0 <= ubound(grd))
    integer,  intent(in) :: iimin   ! min index for ii (>= 0)
    integer,  intent(in) :: iimax   ! max index for ii (<= ubound(grd)+1)
    integer   :: ii  ! result, in [iimin, iimax] 
    !// For ii, grd(ii-1) <= dat < grd(ii), usually. Exceptions are too small/large cases.

    ! Just the initial estimate
    if (dat >= grd(i0 - 1) .and. dat < grd(i0)) then
       ii = i0

       ! Search in the descending order
    else if (dat < grd(i0 - 1)) then
       ii = iimin
       if (dat >= grd(iimin)) then
          do ii = i0 - 1, iimin + 1, -1
             if (dat >= grd(ii - 1)) exit
          end do
       end if

       ! Search in the ascending order
    else
       ii = iimax
       if (dat < grd(iimax - 1)) then ! too large
          do ii = i0 + 1, iimax - 1
             if (dat < grd(ii)) exit
          end do
       end if
    end if

  end function gridIdx_seq0_I


  !+
  ! Search an index in uniform grid
  !-
  function gridIdx_unif0_I(xgrd, xdat, ixmin, ixmax) result(ix)

    real(R_), intent(in) :: xgrd(0:) ! grid values, increasing vector
    real(R_), intent(in) :: xdat     ! target data
    integer,  intent(in) :: ixmin    ! min index
    integer,  intent(in) :: ixmax    ! max index
    integer   :: ix ! result, found index

    ix = int(xdat / (xgrd(1) - xgrd(0))) + 1
    if (ix >= ixmax) then
       ix = ixmax
    else if (xdat < xgrd(ix - 1)) then
       ix = max(ixmin, ix - 1)
    else if (xdat >= xgrd(ix)) then
       ix = ix + 1
    end if

  end function gridIdx_unif0_I


  !+
  ! Grid values in an integer vector
  !    e.g.) J(:) = 2, 4, 6, ..., 12 for (is,iw,n)=(2,2,6)
  !-
  function gridValues_1I(n, is, iw) result(j)

    integer, intent(in) :: n   ! # of grid points to be generated
    integer, intent(in) :: is  ! start point value
    integer, intent(in) :: iw  ! step width
    integer :: j(n)  ! generated grid values
    integer :: i

    do i = 1, n
       j(i) = is + iw * (i - 1)
    end do

  end function gridValues_1I


  !+
  ! Grid values in a real vector
  !    e.g.) rv(:) = (/2, 4, 6, 8, 10/) for (n,rs,rw)=(5,2,2)
  !    e.g.) rv(:) = (/2, 4, 6, 8, 10/) for (n,rs,re)=(5,2,10)
  !-
  function gridValues_1R(n, rs, re, rw) result(rv)

    integer,  intent(in) :: n   ! # of grid points to be generated
    real(R_), intent(in) :: rs  ! start point value
    real(R_), intent(in), optional :: re  ! end point of value
    real(R_), intent(in), optional :: rw  ! step width
    real(R_) :: rv(n)  ! generated grid values
    integer  :: i
    real(R_) :: rw1

    ! Setup
    rw1 = 0.0_R_
    if (present(rw)) rw1 = rw
    if (present(re)) rw1 = (re - rs) / (n - 1) ! possibly overwrite rw

    ! Results
    do i = 1, n
       rv(i) = rs + rw1 * (i - 1)
    end do

  end function gridValues_1R


  !+
  ! Log-linear grid values in a real vector
  !-
  function gridValues_log_1R(n, rs, re) result(rv)

    integer,  intent(in) :: n   ! # of grid points to be generated
    real(R_), intent(in) :: rs  ! start point value
    real(R_), intent(in) :: re  ! end   point value
    real(R_) :: rv(n)  ! generated grid values
    real(R_) :: lnrs, lnre, lnrw
    integer  :: i

    lnrs = log(rs)
    lnre = log(re)
    lnrw = (lnre - lnrs) / real(n - 1, R_)
    do i = 1, n
       rv(i) = exp(lnrs + lnrw * (i - 1))
    end do

  end function gridValues_log_1R


  !+
  ! Average matrix (2D array) data
  !-
  function mat_ave_2R(dat, nx1, ny1, nxw, nyw) result(dat1)

    real(R_), intent(in) :: dat(:,:) ! 2-D original data
    integer,  intent(in) :: nx1, ny1 ! new array sizes
    integer,  intent(in), optional :: nxw, nyw ! # of pixels for averaging
    real(R_)  :: dat1(nx1, ny1) ! result
    integer   :: ix1, iy1, ixs, ixe, iys, iye, nxw1, nyw1
    real(R_)  :: fac

    ! Initialize
    if (present(nxw)) then
       nxw1 = nxw
    else
       nxw1 = size(dat, 1) / nx1
    end if
    if (present(nyw)) then
       nyw1 = nyw
    else
       nyw1 = size(dat, 2) / ny1
    end if
    fac = 1.0_R_ / real(nxw1 * nyw1, R_)

    ! Average
    do iy1 = 1, ny1
       iye = nyw1 * iy1
       iys = iye - nyw1 + 1
       do ix1 = 1, nx1
          ixe = nxw1 * ix1
          ixs = ixe - nxw1 + 1
          dat1(ix1, iy1) = sum(dat(ixs:ixe, iys:iye)) * fac
       end do
    end do

  end function mat_ave_2R


  !+
  ! Matrix operator: Diagonal matrix of a square matrix
  !-
  function mat_diag_2R(w1) result(w2)

    real(R_), intent(in) :: w1(:,:) ! square matrix
    real(R_) :: w2(size(w1,1), size(w1,2))
    integer :: ic

    w2(:,:) = 0.0_R_
    do ic = 1, size(w1,2)
       w2(ic,ic) = w1(ic,ic)
    end do

  end function mat_diag_2R


  !+
  ! Matrix operator: Diagonal matrix of a square matrix in a vector form
  !-
  function mat_diag_1R(w1) result(v2)

    real(R_), intent(in) :: w1(:,:) ! square matrix
    real(R_) :: v2(size(w1,2))
    integer :: ic

    do ic = 1, size(w1,2)
       v2(ic) = w1(ic,ic)
    end do

  end function mat_diag_1R


  !+
  ! Matrix operator: Product of the diagonal elements of a square matrix
  !-
  function mat_diagProd_R(w) result(d)

    real(R_), intent(in) :: w(:,:) ! square matrix
    real(R_) :: d
    integer :: ic

    d = 1.0_R_
    do ic = 1, size(w,2)
       d = d * w(ic,ic)
    end do

  end function mat_diagProd_R

  
  !+
  ! Scale the diagonal elements of a square matrix by a factor
  !-
  subroutine mat_diag_scale(fac, mat)

    real(R_), intent(in) :: fac
    real(R_), intent(inout) :: mat(:,:)
    integer  :: ic
    
    do ic = 1, size(mat, 2)
       mat(ic,ic) = mat(ic,ic) * fac
    end do
       
  end subroutine mat_diag_scale

  
  !+
  ! Scale the nondiagonal elements of a square matrix by a factor
  !-
  subroutine mat_nondiag_scale(fac, mat)

    real(R_), intent(in) :: fac
    real(R_), intent(inout) :: mat(:,:)
    integer  :: ic, nc
    
    nc = size(mat, 2)
    do ic = 1, nc
       if (ic-1 >= 1 ) mat(1:ic-1 , ic) = mat(1:ic-1 , ic) * fac
       if (ic+1 <= nc) mat(ic+1:nc, ic) = mat(ic+1:nc, ic) * fac
    end do
       
  end subroutine mat_nondiag_scale


  !+
  ! Matrix operator: Identity matrix
  !-
  function mat_idt_2R(n) result(w)

    integer, intent(in) :: n
    real(R_) :: w(n,n)
    integer :: i

    w(:,:) = 0.0_R_
    do i = 1, n
       w(i,i) = 1.0_R_
    end do

  end function mat_idt_2R


  !+
  ! Matrix operator: Trace of a square matrix
  !-
  function mat_trace_R(w) result(res)

    real(R_), intent(in) :: w(:,:) ! square matrix
    real(R_) :: res
    integer :: ic

    res = 0.0_R_
    do ic = 1, size(w,2)
       res = res + w(ic,ic)
    end do

  end function mat_trace_R


  !+
  ! Average vector data and reduce data points
  !-
  function vec_ave_1R(vec, nave) result(res) 

    real(R_), intent(in) :: vec(:) ! a vector data
    integer,  intent(in) :: nave   ! # of data points for averaging (should be >= 1)
    real(R_) :: res(size(vec) / nave) ! averaged data with reduced data points
    integer :: nres, ires, is, ie

    nres = size(vec) / nave
    do ires = 1, nres
       is = nave * (ires - 1) + 1
       ie = nave * ires
       res(ires) = sum(vec(is : ie)) / nave
    end do

  end function vec_ave_1R


  !+
  ! Modify a vector by killing too small values
  !  that is, randomly resetting too small values to zero, whereas enlarging some values
  ! Note: All elements of I/O vector and thresholds should be >= 0._R_ 
  ! To save stochastic conservation, the threshold vzero should be pre-determined randomly as
  !    vzero = vmin * RandomNumber
  !-
  subroutine vec_killSmall(vec, vzero, vmin)

    real(R_), intent(inout) :: vec(:) ! I/O vector (should be >= 0)
    real(R_), intent(in)    :: vzero  ! threshold to kill small values
    real(R_), intent(in)    :: vmin   ! min value to survive

    where(vec(:) < vzero)
       vec(:) = 0.0_R_ ! killed
    elsewhere(vec(:) < vmin)
       vec(:) = vmin ! survive
    end where

    !do i = 1, size(vec)
    !   if (vec(i) >= vmin) then ! does nothing
    !   else if (vec(i) < vzero) then
    !      vec(i) = 0.0_R_    ! killed
    !   else          
    !      vec(i) = vmin   ! survive
    !   end if
    !end do
    
  end subroutine vec_killSmall


  !+
  ! Vector operator: Norm of the vector
  !-
  function vec_norm_R(a) result(s) 

    real(R_), intent(in)  :: a(:)    ! vector A
    real(R_) :: s  ! the norm of A
    s = sqrt(sum(a(:)**2))

  end function vec_norm_R


  !+
  ! Outer operation: Compose a matrix from two vectors by the element products
  !-
  function vec_outerSum_2R(a, b) result(w)
    
    real(R_), intent(in) :: a(:), b(:) ! input vectors
    real(R_) :: w(size(a), size(b)) ! output matrix
    w(:,:) = spread(a, dim=2, ncopies=size(b)) + spread(b, dim=1, ncopies=size(a))

  end function vec_outerSum_2R


  !+
  ! Outer operation: Compose a matrix from two vectors by the element products
  !-
  function vec_outerProd_2R(a, b) result(w)
    
    real(R_), intent(in) :: a(:), b(:) ! input vectors
    real(R_) :: w(size(a), size(b)) ! output matrix
    w(:,:) = spread(a, dim=2, ncopies=size(b)) * spread(b, dim=1, ncopies=size(a))

  end function vec_outerProd_2R


  !+
  ! Outer operation: Compose a matrix from two vectors by the element products
  !-
  function vec_outerDiv_2R(a, b) result(w)
    
    real(R_), intent(in) :: a(:), b(:) ! input vectors
    real(R_) :: w(size(a), size(b)) ! output matrix
    w(:,:) = spread(a, dim=2, ncopies=size(b)) / spread(b, dim=1, ncopies=size(a))

  end function vec_outerDiv_2R

end module hparx_vecmat
