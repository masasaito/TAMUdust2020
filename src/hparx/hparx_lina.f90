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
! Library of linear algebra codes
!-
module hparx_lina

  use globals
  use hparx_base,   only : swap, check_idt_I
  use hparx_vecmat, only : mat_diag_1R, mat_diagProd_R, mat_idt_2R, vec_outerProd_2R
  implicit none
  private

  ! Generic interfaces
  interface la_LUDec_solve
     module procedure la_LUDec_solve_m1R, la_LUDec_solve_m2R
  end interface la_LUDec_solve

  ! Public
  public :: la_CholDec
  public :: la_CholDec_inv_2R
  public :: la_CholDec_logDet_R
  public :: la_CholDec_solve_1R
  public :: la_Gauss_inv
  public :: la_Gauss_solve
  public :: la_GaussJ_solve
  public :: la_LUDec
  public :: la_LUDec_solve
  public :: la_LUDec_inv_2R
  public :: la_LUDec_det_R
  public :: la_pivot_row
  public :: la_solveL_1R
  public :: la_solveU_1R
  public :: la_two_inv
  public :: la_two_solve

contains

  !+
  ! Linear Algebra: The Cholesky decomposition of a positive-definite symmetric matrix A
  !-
  subroutine la_CholDec(a, ierr)

    real(R_), intent(inout) :: a(:,:) ! matrix (row, col)
    !// (in)  positive-definite symmetric matrix A(row, col)
    !   (out) Cholesky decomposition (L) of A, where A = L*L^T
    integer,  intent(out) :: ierr     ! error code (0: normal, otherwise: error)
    real(R_) :: tot
    integer  :: n, ic, ir

    ! Initialize
    n = size(a, 1)
    ierr = 0

    ! Lower triangular part with diagonal elements
    do ir = 1, n
       tot = dot_product(a(ir, 1:ir-1), a(ir, 1:ir-1))
       if (a(ir, ir) > tot) then
          a(ir, ir) = sqrt(a(ir, ir) - tot)
          a(ir+1:n, ir) = (a(ir, ir+1:n) - matmul(a(ir+1:n, 1:ir-1), a(ir, 1:ir-1))) / a(ir, ir)
       else
          ierr = 1 ! A is not positive-definite
          return
       end if
    end do

    ! The remaining part
    do ic = 1, n
       a(1:ic-1, ic) = 0.0_R_
    end do

  end subroutine la_CholDec


  !+
  ! Linear Algebra: Determinant of a positive-definite symmetric matrix A
  !  This is just to demonstrate how to compute the determinant by the Cholesky decomposition.
  !-
  function la_CholDec_logDet_R(ell) result(res)

    real(R_), intent(in) :: ell(:,:) ! Cholesky decomposition, L(row,col), of A,
    !// where A = L*L^T is a positive-definite symmetric matrix
    !   Only the diagonal elements are accessed.
    real(R_) :: res ! ln(det), natural logarithm of the determinant

    res = 2.0_R_ * sum(log(mat_diag_1R(ell)))

  end function la_CholDec_logDet_R


  !+
  ! Linear Algebra: Matrix inverse of A using the Cholesky decomposition,
  !  for a positive-definite symmetric matrix A
  !-
  function la_CholDec_inv_2R(ell) result(ainv)

    real(R_), intent(in) :: ell(:,:) ! Cholesky decomposition, L(row,col), of A,
    !// where A = L*L^T is a positive-definite symmetric matrix
    !   Only the lower triangle part of ell and the diagonal elements are accessed.
    real(R_) :: ainv(size(ell,2), size(ell,1)) ! the inverse matrix
    real(R_) :: tot
    integer  :: n, ic, ir

    n = size(ell, 1)

    do ir = 1, n
       do ic = 1, ir
          tot = -dot_product(ell(ir, ic:ir-1), ainv(ic, ic:ir-1))
          if (ir == ic) tot = tot + 1.0_R_
          ainv(ic, ir) = tot / ell(ir, ir)
       end do
    end do

    do ir = n, 1, -1
       do ic = 1, ir
          tot = -dot_product(ell(ir+1:n, ir), ainv(ic, ir+1:n))
          if (ir >= ic) tot = tot + ainv(ic, ir)
          ainv(ir, ic) = tot / ell(ir, ir)
          ainv(ic, ir) = ainv(ir, ic)
       end do
    end do

  end function la_CholDec_inv_2R


  !+
  ! Linear Algebra: Solve a linear matrix equation, A*x = b, where A = L*L^T and
  !  L is the Cholesky decomposition of a positive-definite symmetric matrix A
  !-
  function la_CholDec_solve_1R(ell, b) result(x)

    real(R_), intent(in) :: ell(:,:) ! Cholesky decomposition, L(row,col), of A,
    !// where A = L*L^T is a positive-definite symmetric matrix
    !   Only the lower triangle part of ell and the diagonal elements are accessed.
    real(R_), intent(in) :: b(:) ! a vector in the r.h.s
    real(R_) :: x(size(b)) ! a solution vector
    integer  :: n, ir

    ! Solve L*x' = b
    x(:) = la_solveL_1R(ell, b)

    ! Solve L^T*x = x'
    !x(:) = la_solveU_1R(transpose(ell), x)
    n = size(ell, 1)
    do ir = n, 1, -1 ! bottom to top
       x(ir) = (x(ir) - dot_product(ell(ir+1:n, ir), x(ir+1:n))) / ell(ir, ir)
    end do

  end function la_CholDec_solve_1R


  !+
  ! Linear Algebra: Inverse matrix A^-1, for a square matrix A by the Gaussian elimination
  !-
  subroutine la_Gauss_inv(a, ainv, ierr)

    real(R_), intent(inout) :: a(:,:) ! square matrix (row, col)
    !// (in)  : input matrix A
    !// (out) : upper triangular matrix transformed from A
    real(R_), intent(out) :: ainv(:,:) ! inverse of matrix A
    integer,  intent(out) :: ierr      ! error code (0: normal, otherwise: error)

    ainv(:,:) = mat_idt_2R(size(a,1))
    call la_Gauss_solve(a, ainv, ierr)

  end subroutine la_Gauss_inv


  !+
  ! Linear Algebra: Solve linear matrix equations, A*x=b, by the Gaussian elimination
  !  with backsubstitution and partial pivoting
  !-
  subroutine la_Gauss_solve(a, b, ierr)

    real(R_), intent(inout) :: a(:,:) ! square matrix (row, col)
    !// (in)  : input matrix A
    !// (out) : upper triangular matrix transformed from A
    real(R_), intent(inout) :: b(:,:) ! matrix with vectors in columns
    !// (in)  : input vectors b's in the r.h.s
    !   (out) : solution vectors x's
    !   where A(:,:)*x(:,i) = b(:,i), for all i
    integer,  intent(out) :: ierr     ! error code (0: normal, otherwise: error)
    integer  :: n, ic, ir, iv, nv
    real(R_) :: apiv

    ! Initialize
    n  = size(a, 1)
    nv = size(b, 2)
    ierr = 0

    ! Forward elimination with partial pivoting
    do ic = 1, n
       call la_pivot_row(a, ic, ir, ierr)
       if (ierr /= 0) return
       if (ir /= ic) then
          call swap(a(ic, ic:n), a(ir, ic:n))
          call swap(b(ic, :), b(ir, :))
       end if
       apiv = 1.0_R_ / a(ic, ic) ! adjoint of the pivot
       a(ic, ic) = 1.0_R_
       a(ic, ic+1:n) = a(ic, ic+1:n) * apiv
       b(ic, :)      = b(ic, :)      * apiv
       do ir = ic + 1, n
          a(ir, ic+1:n) = a(ir, ic+1:n) - a(ir, ic) * a(ic, ic+1:n)
          b(ir, :)      = b(ir, :)      - a(ir, ic) * b(ic, :)
          a(ir, ic) = 0.0_R_
       end do
    end do

    ! Backward substitution
    do iv = 1, nv
       do ir = n, 1, -1
          b(ir, iv) = b(ir, iv) - dot_product(a(ir, ir+1:n), b(ir+1:n, iv))
       end do
    end do

  end subroutine la_Gauss_solve


  !+
  ! Linear Algebra: Solve linear matrix equations, A*x=b, by the Gaussian-Jordan elimination
  !  with partial pivoting
  ! Note : The Gaussian elimination with backsubstitution outperforms.
  !-
  subroutine la_GaussJ_solve(a, b, ierr)

    real(R_), intent(inout) :: a(:,:) ! square matrix (row, col)
    !// (in)  : input matrix A
    !// (out) : identity matrix transformed from A
    real(R_), intent(inout) :: b(:,:) ! matrix with vectors in columns
    !// (in)  : input vectors b's in the r.h.s
    !   (out) : solution vectors x's
    !   where A(:,:)*x(:,i) = b(:,i), for all i
    integer,  intent(out) :: ierr     ! error code (0: normal, otherwise: error)
    integer  :: n, ic, ir
    real(R_) :: apiv

    ! Initialize
    n = size(a, 1)
    ierr = 0

    ! Elimination with partial pivoting
    do ic = 1, n
       call la_pivot_row(a, ic, ir, ierr)
       if (ierr /= 0) return
       if (ir /= ic) then
          call swap(a(ic, ic:n), a(ir, ic:n))
          call swap(b(ic, :), b(ir, :))
       end if
       apiv = 1.0_R_ / a(ic, ic) ! adjoint of the pivot
       a(ic, ic) = 1.0_R_
       a(ic, ic+1:n) = a(ic, ic+1:n) * apiv
       b(ic, :)      = b(ic, :)      * apiv
       do ir = 1, n
          if (ir /= ic) then
             a(ir, ic+1:n) = a(ir, ic+1:n) - a(ir, ic) * a(ic, ic+1:n)
             b(ir, :)      = b(ir, :)      - a(ir, ic) * b(ic, :)
             a(ir, ic) = 0.0_R_
          end if
       end do
    end do

  end subroutine la_GaussJ_solve


  !+
  ! Linear Algebra: LU decomposition of a matrix
  !  The Crout's algorithm (with partial pivoting)
  !
  ! Reference: NR 3rd Ed, p.48-51
  !-
  subroutine la_LUDec(a, iir, odd, ierr)

    real(R_), intent(inout) :: a(:,:) ! square matrix (row, col)
    !// (in)  : input matrix
    !   (out) : LU decomposition (non-diagonal part of L and full U)
    integer,  intent(out) :: iir(:) ! row index table (row-interchanging information)
    logical,  intent(out) :: odd    ! true if # of interchanges is odd
    integer,  intent(out) :: ierr   ! error code (0: normal, otherwise: error)
    real(R_), parameter :: SML = 1.0e-20_R_
    real(R_) :: fac(size(a,1))
    integer  :: n, ic, ir

    ! Initialize
    n = check_idt_I(size(a,1), size(a,2), size(iir), title='la_LUDec: ')
    ierr = 0
    odd = .false.

    ! Get the scaling factors
    fac(:) = maxval(abs(a(:,:)), dim=2)
    if (any(fac(:) == 0.0_R_)) then
       ierr = 1 ! singular matrix with a row of zeros
       return
    end if
    fac(:) = 1.0_R_ / fac(:)

    ! Crout's algorithm
    do ic = 1, n
       ir = ic - 1 + maxloc(fac(ic:n) * abs(a(ic:n,ic)), dim=1) ! index for the largest pivot element
       if (ir >= ic .and. ir <= n) then
          iir(ic) = ir ! save the index for output
          if (ic /= ir) then ! interchange rows
             call swap(a(ir,:), a(ic,:))
             fac(ir) = fac(ic)
             if (odd) then
                odd = .false.
             else
                odd = .true.
             end if
          end if
          if (a(ic,ic) == 0.0_R_) a(ic,ic) = SML
          if (ic < n) then ! this clause is nessesary to avoid SIGSEG.
             a(ic+1:n, ic) = a(ic+1:n, ic) / a(ic,ic)
             a(ic+1:n, ic+1:n) = a(ic+1:n, ic+1:n) - vec_outerProd_2R(a(ic+1:n, ic), a(ic, ic+1:n))
          end if
          !// Reduce remaining submatrix: Outer product Gaussian elimination by Golub and Van Loan
       else
          print *, 'la_LUDec: Error. The matrix has NaN or Inf?'
          print *, n, ir, ic, maxloc(fac(ic:n) * abs(a(ic:n,ic)), dim=1)
          print *, a(:,:)
          ierr = 2 ! matrix with infinite values or NAN
          return
       endif
    end do

  end subroutine la_LUDec


  !+
  ! Linear Algebra: Solve the coupled lienar equations, A*x=b, by using the LU decomposition
  !  Assuming that the partial pivoting is operated.
  !
  ! Reference: NR 3rd Ed, p.48-49
  !-
  subroutine la_LUDec_solve_m1R(alu, iir, b)

    real(R_), intent(in) :: alu(:,:) ! square matrix (row, col)
    !// LU decomposition of A, with the non-diagonal part of L and full U
    integer,  intent(in) :: iir(:) ! row index table (row-interchanging information)
    real(R_), intent(inout) :: b(:) ! vector
    !// (in)  the vector b
    !   (out) the vector x, the solution
    real(R_) :: wrk
    integer :: n, ir, irb, irp

    ! Forward substitution
    n = check_idt_I(size(alu,1), size(alu,2), size(b), size(iir), title='la_LUDec_solve_m1R: ')
    irb = 0
    do ir = 1, n
       irp = iir(ir) ! to unscramble the permutation
       wrk = b(irp)
       b(irp) = b(ir)
       if (irb /= 0) then
          wrk = wrk - dot_product(alu(ir, irb:ir-1), b(irb:ir-1))
       else if (wrk /= 0.0_R_) then ! nonzero element was detected
          irb = ir ! index of the first nonzero element of b
       end if
       b(ir) = wrk
    end do

    ! Backsubstitution
    b(n) = b(n) / alu(n,n)
    do ir = n-1, 1, -1
       b(ir) = (b(ir) - dot_product(alu(ir, ir+1:n), b(ir+1:n))) / alu(ir,ir)
    end do

  end subroutine la_LUDec_solve_m1R


  !+
  ! Linear Algebra: Solve the coupled lienar equations, A*X=B, by using the LU decomposition
  !  Assuming that the partial pivoting is operated.
  !
  ! Reference: NR 3rd Ed, p.48-49
  !-
  subroutine la_LUDec_solve_m2R(alu, iir, b)

    real(R_), intent(in) :: alu(:,:) ! N*N square matrix (row, col)
    !// LU decomposition of A, with the non-diagonal part of L and full U
    integer,  intent(in) :: iir(:) ! row index table (row-interchanging information)
    real(R_), intent(inout) :: b(:,:) ! N*M matrix (row, col)
    !// (in)  the matrix B
    !   (out) the matrix X, the solution
    integer :: ic

    do ic = 1, size(b,2)
       call la_LUDec_solve_m1R(alu, iir, b(:,ic))
    end do

  end subroutine la_LUDec_solve_m2R


  !+
  ! Linear Algebra: Solve the coupled lienar equation, A*X=B, by using the LU decomposition
  !  Assuming that the partial pivoting is operated.
  !
  ! Reference: NR 3rd Ed, p.54
  !-
  function la_LUDec_inv_2R(alu, iir) result(ainv)

    real(R_), intent(in) :: alu(:,:) ! square matrix (row, col)
    !// LU decomposition of A, with the non-diagonal part of L and full U
    integer,  intent(in) :: iir(:) ! row index table (row-interchanging information)
    real(R_) :: ainv(size(alu,1), size(alu,1)) ! the inverse matrix (row, col)

    ainv(:,:) = mat_idt_2R(size(alu,1))
    call la_LUDec_solve_m2R(alu, iir, ainv)

  end function la_LUDec_inv_2R


  !+
  ! Linear Algebra: Determinant of the matrix A, using the LU decomposition
  !
  ! Reference: NR 3rd Ed, p.55
  ! Note: Overflow or underflow may occur easily when N is large.
  !    A modified algorithm should be used in such cases.
  !-
  function la_LUDec_det_R(alu, odd) result(res)

    real(R_), intent(in) :: alu(:,:) ! square matrix (row, col)
    !// LU decomposition of A, with the non-diagonal part of L and full U
    logical,  intent(in) :: odd  ! true if # of interchanges is odd
    real(R_) :: res ! determinant of A

    res = mat_diagProd_R(alu)
    if (odd) res = -res
    !// This is simply to demonstrate how to get the determinant from the LU decomposition

  end function la_LUDec_det_R


  !+
  ! Linear Algebra: Partial pivoting, interchanging rows of a matrix
  !-
  subroutine la_pivot_row(a, ip, ir, ierr)

    real(R_), intent(in)  :: a(:,:) ! square matrix (row, col)
    integer,  intent(in)  :: ip     ! present pivot index
    integer,  intent(out) :: ir     ! replaced row index
    integer,  intent(out) :: ierr   ! error code (0=normal)
    integer :: n

    n = size(a, 1)
    ierr = 0
    ir = ip - 1 + maxloc(abs(a(ip:n, ip)), 1) ! find the maximum
    if (a(ir, ip) == 0.0_R_) then ! A is singular
       ierr = 1
       return
    end if

  end subroutine la_pivot_row


  !+
  ! Linear Algebra: Solve a linear matrix equation, L*x=b, where L is a lower triangle matrix
  !-
  function la_solveL_1R(ell, b) result(x)

    real(R_), intent(in) :: ell(:,:) ! lower triangle matrix L(row,col)
    !// Only the lower triangle part of ell and the diagonal elements are accessed.
    !   The diagonal elements should be non-zero.
    real(R_), intent(in) :: b(:) ! a vector in the r.h.s
    real(R_) :: x(size(b)) ! a solution vector
    integer  :: n, ir

    n = size(ell, 1)
    do ir = 1, n ! top to bottom
       x(ir) = (b(ir) - dot_product(ell(ir, 1:ir-1), x(1:ir-1))) / ell(ir, ir)
    end do

  end function la_solveL_1R


  !+
  ! Linear Algebra: Solve a linear matrix equation, U*x=b, where U is a upper triangle matrix
  !-
  function la_solveU_1R(u, b) result(x)

    real(R_), intent(in) :: u(:,:) ! upper triangle matrix U(row,col)
    !// Only the upper triangle part of ell and the diagonal elements are accessed.
    !   The diagonal elements should be non-zero.
    real(R_), intent(in) :: b(:) ! a vector in the r.h.s
    real(R_) :: x(size(b)) ! a solution vector
    integer  :: n, ir

    n = size(u, 1)
    do ir = n, 1, -1 ! bottom to top
       x(ir) = (b(ir) - dot_product(u(ir, ir+1:n), x(ir+1:n))) / u(ir, ir)
    end do

  end function la_solveU_1R


  !+
  ! Linear Algebra: Solve linear matrix equations, A*x=b, for a 2x2 matrix A
  !-
  subroutine la_two_inv(a, ainv, ierr)

    real(R_), intent(in)  :: a(:,:)    ! square 2x2 matrix A(row, col)
    real(R_), intent(out) :: ainv(:,:) ! inverse of matrix A
    integer,  intent(out) :: ierr      ! error code (0: normal, otherwise: error)
    real(R_) :: b(2,2)

    b(1:2, 1) = (/ 1.0_R_, 0.0_R_ /)
    b(1:2, 2) = (/ 0.0_R_, 1.0_R_ /)
    call la_two_solve(a, b, ainv, ierr)

  end subroutine la_two_inv


  !+
  ! Linear Algebra: Solve linear matrix equations, A*x=b, for a 2x2 matrix A
  !-
  subroutine la_two_solve(a, b, x, ierr)

    real(R_), intent(in)  :: a(:,:) ! square 2x2 matrix A(row, col)
    real(R_), intent(in)  :: b(:,:) ! matrix with input vectors b's in columns
    real(R_), intent(out) :: x(:,:) ! matrix with solution vectors x's in columns
    !// where A(:,:)*x(:,i) = b(:,i), for all i
    integer,  intent(out) :: ierr   ! error code (0: normal, otherwise: error)
    real(R_) :: aa11, d, dinv

    ! Initialize
    ierr = 0
    d = a(1, 1) * a(2, 2) - a(2, 1) * a(1, 2)
    if (d == 0.0_R_) then ! A is singular
       ierr = 1
       return
    end if
    dinv = 1.0_R_ / d

    ! No pivoting
    if (abs(a(1,1)) >= abs(a(2,1))) then
       if (a(1, 1) == 0.0_R_) then ! A is singular
          ierr = 1
          return
       end if
       aa11 = 1.0_R_ / a(1, 1) ! adjoint of the pivot
       x(2, :) = (a(1, 1) * b(2, :) - a(2, 1) * b(1, :)) * dinv
       x(1, :) = (          b(1, :) - a(1, 2) * x(2, :)) * aa11

       ! With pivoting
    else
       if (a(2, 1) == 0.0_R_) then ! A is singular
          ierr = 1
          return
       end if
       aa11 = 1.0_R_ / a(2, 1) ! adjoint of the pivot
       x(2, :) = (a(1, 1) * b(2, :) - a(2, 1) * b(1, :)) * dinv
       x(1, :) = (          b(2, :) - a(2, 2) * x(2, :)) * aa11
    end if

  end subroutine la_two_solve

end module hparx_lina
