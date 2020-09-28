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
! Library of basic utilities
!-
module hparx_base

  use globals, only : R_, RD_, R4_, R8_, I4_, I8_, IMIN_, IMAX_, RLRG_
  implicit none
  private

  ! Generic interfaces
  interface swap
     module procedure swap_mI, swap_mR, swap_m1R, swap_m2R
  end interface swap

  ! Public
  public :: abstract_1R
  public :: abstract_2R
  public :: bswap_m1I
  public :: bswap_m1R4
  public :: bswap_m1A
  public :: calen_addDays
  public :: calen_JD_R
  public :: calen_MJD_R
  public :: calen_month_A3
  public :: calen_month_I
  public :: calen_jDay_I
  public :: calen_lagDays_I
  public :: calen_leapYear_I
  public :: check_args
  public :: check_idt_I
  public :: check_i1I
  public :: check_i1R
  public :: check_iI
  public :: check_iR
  public :: check_valid
  public :: convLtl_AN_to_I
  public :: convLtl_AN_to_1I
  public :: cyclicAdd_I
  public :: err_setUnit
  public :: err_setAct
  public :: err_issue
  public :: err_open
  public :: err_read
  public :: err_write
  public :: fileName_AN
  public :: filePath_AN
  public :: freeUnit_I
  public :: num2str_AN
  public :: nonZero_R
  public :: numDigits_I
  public :: read1word
  public :: shuffle_m1I
  public :: shuffle_m1A
  public :: sort_heap
  public :: sort_selec
  public :: str2num_I
  public :: str2num_skip_I
  public :: str_afterStr_AN
  public :: str_shift_AN
  public :: str_split
  public :: swap
  public :: writeStr

  ! Private variables
  integer, save :: Err_iue = 6  ! unit index for standard error output
  integer, save :: Err_mact = 1 ! action to do when an error is detected
  ! - mact =  1, the execution will stop after printing some error messages.
  ! - mact =  0, just prints warning messages, and continue the execution.
  ! - mact = -1, does nothing even if some fatal error might is detected.

contains

  !+
  ! Abstract (extract) 1-D real array
  ! ex) If (is, iw, n) = (3, 2, 4), then retrieved points will be (3, 5, 7, 9)
  !-
  function abstract_1R(dat, n1, is, iw) result(dat1)

    real(R_), intent(in) :: dat(:) ! original data
    integer,  intent(in) :: n1     ! # of points retrieved
    integer,  intent(in), optional :: is ! start point
    integer,  intent(in), optional :: iw ! spacing width in points
    real(R_) :: dat1(n1) ! result
    integer :: is1, ie1, iw1, n2

    is1 = 1
    iw1 = 1
    if (present(is)) is1 = is ! start
    if (present(iw)) iw1 = iw ! step width

    ie1 = is1 + iw1 * (n1 - 1)
    n2 = n1
    if (ie1 > size(dat)) then ! exceed the upper bound
       n2 = (size(dat) - is1) / iw1 + 1
       ie1 = is1 + iw1 * (n2 - 1)
    end if

    dat1(1:n2) = dat(is1:ie1:iw1) ! conversion

  end function abstract_1R


  !+
  ! Abstract (extract) 2-D real array
  ! ex) If (ixs, ixw, nx) = (3, 2, 4), then retrieved X points will be (3, 5, 7, 9)
  !-
  function abstract_2R(dat, nx1, ny1, ixs, ixw, iys, iyw) result(dat1)

    real(R_), intent(in) :: dat(:,:)  ! original data
    integer,  intent(in) :: nx1, ny1  ! # of X/Y points retrieved
    integer,  intent(in), optional :: ixs ! X start point
    integer,  intent(in), optional :: ixw ! X spacing width in points
    integer,  intent(in), optional :: iys ! Y start point
    integer,  intent(in), optional :: iyw ! Y spacing width in points
    real(R_) :: dat1(nx1, ny1) ! result
    integer :: ixs1, ixe1, ixw1, iys1, iye1, iyw1, nx2, ny2

    ! Setup for X
    ixs1 = 1
    ixw1 = 1
    if (present(ixs)) ixs1 = ixs ! start
    if (present(ixw)) ixw1 = ixw ! step width
    ixe1 = ixs1 + ixw1 * (nx1 - 1)
    nx2 = nx1
    if (ixe1 > size(dat, 1)) then ! exceed the upper bound
       nx2 = (size(dat, 1) - ixs1) / ixw1 + 1
       ixe1 = ixs1 + ixw1 * (nx2 - 1)
    end if

    ! Setup for Y
    iys1 = 1
    iyw1 = 1
    if (present(iys)) iys1 = iys ! start
    if (present(iyw)) iyw1 = iyw ! step width
    iye1 = iys1 + iyw1 * (ny1 - 1)
    ny2 = ny1
    if (iye1 > size(dat, 2)) then ! exceed the upper bound
       ny2 = (size(dat, 2) - iys1) / iyw1 + 1
       iye1 = iys1 + iyw1 * (ny2 - 1)
    end if

    ! Conversion
    dat1(1:nx2, 1:ny2) = dat(ixs1:ixe1:ixw1, iys1:iye1:iyw1)

  end function abstract_2R


  !+
  ! Byte swap for integer 1-D array
  !-
  subroutine bswap_m1I(jj, nx, nb)

    integer, intent(inout) :: jj(:)      ! data vector
    integer, intent(in), optional :: nx  ! # of data processed
    integer, intent(in), optional :: nb  ! # of bytes for swapping [2/4/8] (default = 4)
    integer, parameter :: NWRKMAX = 100000  ! max size of work vector
    integer :: nx1, nb1, na, nc, nwrk, ijob, njob, ixs, ixe
    character(1), save :: at(1) = ' '
    character(1), allocatable :: aa(:), a(:)

    ! Optional arguments
    nx1 = size(jj)
    nb1 = 4
    if (present(nx)) nx1 = nx
    if (present(nb)) nb1 = nb
    if (nb1 /= 2 .and. nb1 /= 4 .and. nb1 /= 8) return ! invalid nb

    ! Setup
    nwrk = min(NWRKMAX, nx1)       ! size of work vector
    nwrk = (nwrk / nb1 + 1) * nb1  ! reset to multiple of nb1
    njob = (nx1 - 1) / nwrk + 1    ! # of jobs
    nc = size(transfer(jj(1), at)) ! # of characters for single data
    na = nc * nwrk                 ! total # of characters for single job
    allocate (aa(na), a(na / nb1))

    ! Do jobs
    do ijob = 1, njob
       ixs = nwrk * (ijob - 1) + 1
       ixe = min(nx1, nwrk * ijob)
       na = nc * (ixe - ixs + 1)
       aa(1:na) = transfer(jj(ixs:ixe), aa)
       call bswap_m1A(aa, a, nb1)
       jj(ixs:ixe) = transfer(aa(1:na), jj)
    end do
    deallocate (aa, a)

  end subroutine bswap_m1I


  !+
  ! Byte swap for real(R4_) 1-D array
  !-
  subroutine bswap_m1R4(rr, nx, nb)

    real(R4_), intent(inout) :: rr(:)      ! data vector
    integer,   intent(in), optional :: nx  ! # of data processed
    integer,   intent(in), optional :: nb  ! # of bytes for swapping [2/4/8] (default = 4)
    integer, parameter :: NWRKMAX = 100000  ! max size of work vector
    integer :: nx1, nb1, na, nc, nwrk, ijob, njob, ixs, ixe
    character(1), save :: at(1) = ' '
    character(1), allocatable :: aa(:), a(:)

    ! Optional arguments
    nx1 = size(rr)
    nb1 = 4
    if (present(nx)) nx1 = nx
    if (present(nb)) nb1 = nb
    if (nb1 /= 2 .and. nb1 /= 4 .and. nb1 /= 8) return ! invalid nb

    ! Setup
    nwrk = min(NWRKMAX, nx1)       ! size of work vector
    nwrk = (nwrk / nb1 + 1) * nb1  ! reset to multiple of nb1
    njob = (nx1 - 1) / nwrk + 1    ! # of jobs
    nc = size(transfer(rr(1), at)) ! # of characters for single data
    na = nc * nwrk                 ! total # of characters for single job
    allocate (aa(na), a(na / nb1))

    ! Do jobs
    do ijob = 1, njob
       ixs = nwrk * (ijob - 1) + 1
       ixe = min(nx1, nwrk * ijob)
       na = nc * (ixe - ixs + 1)
       aa(1:na) = transfer(rr(ixs:ixe), aa)
       call bswap_m1A(aa, a, nb1)
       rr(ixs:ixe) = transfer(aa(1:na), rr)
    end do
    deallocate (aa, a)

  end subroutine bswap_m1R4


  !+
  ! Byte swap for a character(1) vector
  !-
  subroutine bswap_m1A(adat, awrk, nb, nwrk)

    character(1), intent(inout) :: adat(:) ! target data vector
    character(1), intent(inout) :: awrk(:) ! work vector
    integer,      intent(in)    :: nb      ! # of bytes for swapping [2/4/8]
    integer,      intent(in), optional :: nwrk  ! # of swapping works
    integer :: nw, na

    ! Sizes
    nw = size(awrk)
    if (present(nwrk)) nw = nwrk
    na = nb * nw

    ! Swap
    if (nb == 2) then         ! 2 bytes
       awrk(1:nw)      = adat(1:na-1:nb)
       adat(1:na-1:nb) = adat(2:na  :nb)
       adat(2:na  :nb) = awrk(1:nw)
    else if (nb == 4) then    ! 4 bytes
       awrk(1:nw)      = adat(1:na-3:nb)
       adat(1:na-3:nb) = adat(4:na  :nb)
       adat(4:na  :nb) = awrk(1:nw)
       awrk(1:nw)      = adat(2:na-2:nb)
       adat(2:na-2:nb) = adat(3:na-1:nb)
       adat(3:na-1:nb) = awrk(1:nw)
    else if (nb == 8) then    ! 8 bytes
       awrk(1:nw)      = adat(1:na-7:nb)
       adat(1:na-7:nb) = adat(8:na  :nb)
       adat(8:na  :nb) = awrk(1:nw)
       awrk(1:nw)      = adat(2:na-6:nb)
       adat(2:na-6:nb) = adat(7:na-1:nb)
       adat(7:na-1:nb) = awrk(1:nw)
       awrk(1:nw)      = adat(3:na-5:nb)
       adat(3:na-5:nb) = adat(6:na-2:nb)
       adat(6:na-2:nb) = awrk(1:nw)
       awrk(1:nw)      = adat(4:na-4:nb)
       adat(4:na-4:nb) = adat(5:na-3:nb)
       adat(5:na-3:nb) = awrk(1:nw)
    end if

  end subroutine bswap_m1A


  !+
  ! Add days to a date
  !-
  subroutine calen_addDays(iyyyy, imm, idd, ndays)

    integer,   intent(inout) :: iyyyy ! year
    integer,   intent(inout) :: imm   ! month
    integer,   intent(inout) :: idd   ! day
    integer,   intent(in) :: ndays    ! # of days (larger than 0)
    integer   :: ijul, ileap, njul
    integer, save :: ndd(0:12) = (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365/)

    ! Days from 00Z Jan 1st
    ijul = ndd(imm - 1) + idd
    if (imm >= 3) ijul = ijul + calen_leapYear_I(iyyyy)

    ! Add days
    ijul = ijul + ndays
    do
       njul = 365 + calen_leapYear_I(iyyyy)
       if (ijul <= njul) exit
       iyyyy = iyyyy + 1
       ijul = ijul - njul
    end do

    ! Retrieve month & day
    if (ijul <= ndd(1)) then
       imm = 1
       idd = ijul
    else
       ileap = calen_leapYear_I(iyyyy)
       do imm = 2, 12
          if (ijul <= ndd(imm) + ileap) exit
       end do
       idd = ijul - ndd(imm - 1)
    end if

  end subroutine calen_addDays


  !+
  ! Julian date (JD)
  !-
  function calen_JD_R(iy, im, day) result(res)

    integer, intent(in) :: iy, im ! year, month
    real(R_), intent(in) :: day ! date (12.5 for 12:00 UTC of day 12)
    real(R_) :: res ! JD
    res = calen_MJD_R(iy, im, day) + 2400000.5_R_

  end function calen_JD_R


  !+
  ! Modified Julian date (MJD)
  !-
  function calen_MJD_R(iy, im, day) result(res)

    integer,  intent(in) :: iy, im ! year, month
    real(R_), intent(in) :: day ! date (12.5 for 12:00 UTC of day 12)
    real(R_) :: res ! MJD
    !// The Gregorian calendar is assumed after Sep. 3, 1752, which is common in UNIX systems.
    !   Note the Gregorian calendar was proposed on Oct. 15, 1582 and the United Kingdom
    !   adopted the Gregorian in 1752.
    integer :: iyy, imm

    ! Initialize
    iyy = iy
    imm = im
    if (imm <= 2) then
       iyy = iyy - 1
       imm = imm + 12
    end if
    
    ! MJD
    if (iy >= 1753 .or. (iy == 1752 .and. im >= 10) .or. &
         (iy == 1752 .and. im == 9 .and. day >= 3.0_R_)) then ! the Gregorian calendar
       res = int(365.25_R_ * iyy) + int(iyy / 400.0_R_) - int(iyy / 100.0_R_) &
            + int(30.59_R_ * (imm - 2)) + day - 678912.0_R_
    else if (iyy >= 0) then ! AC
       res = int(365.25_R_ * iyy) + int(30.59_R_ * (imm - 2)) + day - 678914.0_R_
    else ! BC
       res = int(365.25_R_ * iyy) + int(30.59_R_ * (imm - 2)) + day - 678915.0_R_
    end if

  end function calen_MJD_R


  !+
  ! Month character string (e.g., 'DEC') from a month number
  !-
  function calen_month_A3(imt) result(smt)

    integer, intent(in) :: imt ! month numeric value (1-12)
    character(3)        :: smt ! month character string ('DEC' for December)
    character(len=36), save :: mttab = 'JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC'
    smt = mttab(3*imt-2 : 3*imt)

  end function calen_month_A3


  !+
  ! Month number from a character string
  !-
  function calen_month_I(smt) result(imt)

    character(*), intent(in) :: smt ! month character string ('DEC' for December)
    integer                  :: imt ! month numeric value (1-12)
    character(len=36), save :: mttab = 'JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC'
    imt = (index(mttab, smt(1:3)) + 2) / 3

  end function calen_month_I


  !+
  ! Day of Year (DOY), days (1-366) counted from the beginning of year.
  !-
  function calen_jDay_I(iyr, imon, iday) result(jday)

    integer, intent(in) :: iyr, imon, iday ! year, month, day
    integer :: jday ! DOY
    integer, save :: maxday(12) = (/0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334/)
    jday = maxday(imon) + iday
    if (imon >= 3) jday = jday + calen_leapYear_I(iyr)

  end function calen_jDay_I


  !+
  ! Lag of days between given two dates
  !-
  function calen_lagDays_I(iy0, im0, id0, iy1, im1, id1) result(jd)

    integer, intent(in) :: iy0, im0, id0
    integer, intent(in) :: iy1, im1, id1 ! iy1 should be >= iy0
    integer :: jd, iytmp

    jd = 0
    if (iy1 == iy0) then
       jd = calen_jDay_I(iy1, im1, id1) - calen_jDay_I(iy0, im0, id0)
    else if (iy1 > iy0) then
       jd = calen_jDay_I(iy0, 12, 31) - calen_jDay_I(iy0, im0, id0)
       do iytmp = iy0 + 1, iy1 - 1
          jd = jd + calen_jDay_I(iytmp, 12, 31)
       end do
       jd = jd + calen_jDay_I(iy1, im1, id1)
    else
       jd = -2000000000
    endif

  end function calen_lagDays_I


  !+
  ! Leap year or not
  !-
  function calen_leapYear_I(iyyyy) result(i)

    integer,   intent(in) :: iyyyy ! year
    integer   :: i  ! result: 1 for leap year, 0 for non-leap year

    if (mod(iyyyy, 400) == 0) then
       i = 1
    else if (mod(iyyyy, 100) == 0) then
       i = 0
    else if (mod(iyyyy, 4) == 0) then
       i = 1
    else
       i = 0
    end if

  end function calen_leapYear_I


  !+
  ! Check arguments and possibly dies with an error message if any logical is false
  !-
  subroutine check_args(n1, n2, n3, n4, n5, title)

    logical, intent(in) :: n1 ! first argument
    logical, optional, intent(in) :: n2, n3, n4, n5 ! optional arguments
    character(len=*), optional, intent(in) :: title ! title error message
    character(len=21) :: tag = 'check_args failed :'

    if (.not. n1) call err_issue(1, title//tag)
    if (present(n2)) then
       if (.not. n2) call err_issue(1, title//tag)
       if (present(n3)) then
          if (.not. n3) call err_issue(1, title//tag)
          if (present(n4)) then
             if (.not. n4) call err_issue(1, title//tag)
             if (present(n5)) then
                if (.not. n5) call err_issue(1, title//tag)
             end if
          end if
       end if
    end if

  end subroutine check_args


  !+
  ! Check the identity of multiple integer arguments and possibly dies with an error message
  !-
  function check_idt_I(n1, n2, n3, n4, n5, title) result(n)

    integer, intent(in) :: n1, n2 ! first and second arguments
    integer, optional, intent(in) :: n3, n4, n5 ! optional arguments
    character(len=*), optional, intent(in) :: title ! title error message
    character(len=20) :: tag = 'check_idt_I failed :'
    integer :: n

    if (n1 .ne. n2) call err_issue(1, title//tag)
    if (present(n3)) then
       if (n1 .ne. n3) call err_issue(1, title//tag)
       if (present(n4)) then
          if (n1 .ne. n4) call err_issue(1, title//tag)
          if (present(n5)) then
             if (n1 .ne. n5) call err_issue(1, title//tag)
          end if
       end if
    end if
    n = n1

  end function check_idt_I


  !+
  ! Check whether integer vector values is in a valide range
  !-
  subroutine check_i1I(vname, jval, imin, imax)

    character(*), intent(in) :: vname   ! variable name (possibly with arbitrary comments)
    integer,      intent(in) :: jval(:) ! value to be checked
    integer,      intent(in), optional :: imin ! min limit
    integer,      intent(in), optional :: imax ! max limit
    ! - Note: Should be imin <= jval <= imax, in a normal case.
    character(256) :: msg
    integer :: i

    do i = 1, size(jval)
       if (present(imin)) then
          if (jval(i) < imin) then
             write (msg, *) trim(vname)//' = ', jval(i), ': Too small. Should be >= ', imin
             call err_issue(-i, msg)
          end if
       end if

       if (present(imax)) then
          if (jval(i) > imax) then
             write (msg, *) trim(vname)//' = ', jval(i), ': Too large. Should be <= ', imax
             call err_issue(i, msg)
          end if
       end if
    end do

  end subroutine check_i1I


  !+
  ! Check whether real vector values is in a valide range
  !-
  subroutine check_i1R(vname, rval, rmin, rmax)

    character(*), intent(in) :: vname   ! variable name (possibly with arbitrary comments)
    real(R_),     intent(in) :: rval(:) ! value to be checked
    real(R_),     intent(in), optional :: rmin ! min limit
    real(R_),     intent(in), optional :: rmax ! max limit
    ! - Note: Should be rmin <= rval <= rmax, in a normal case.
    character(256) :: msg
    integer :: i

    do i = 1, size(rval)
       if (present(rmin)) then
          if (rval(i) < rmin) then
             write (msg, *) trim(vname)//' = ', rval(i), ': Too small. Should be >= ', rmin
             call err_issue(-i, msg)
          end if
       end if

       if (present(rmax)) then
          if (rval(i) > rmax) then
             write (msg, *) trim(vname)//' = ', rval(i), ': Too large. Should be <= ', rmax
             call err_issue(i, msg)
          end if
       end if
    end do

  end subroutine check_i1R


  !+
  ! Check whether an integer scalar value is in a valid range
  !-
  subroutine check_iI(vname, ival, imin, imax)

    character(*), intent(in) :: vname ! variable name (possibly with arbitrary comments)
    integer,      intent(in) :: ival  ! value to be checked
    integer,      intent(in), optional :: imin ! min limit
    integer,      intent(in), optional :: imax ! max limit
    ! - Note: Should be imin <= ival <= imax, in a normal case.
    character(256) :: msg

    if (present(imin)) then
       if (ival < imin) then
          write (msg, *) trim(vname)//' = ', ival, ': Too small. Should be >= ', imin
          call err_issue(-1, msg)
       end if
    end if

    if (present(imax)) then
       if (ival > imax) then
          write (msg, *) trim(vname)//' = ', ival, ': Too large. Should be <= ', imax
          call err_issue(1, msg)
       end if
    end if

  end subroutine check_iI


  !+
  ! Check whether a real scalar value is in a valide range
  !-
  subroutine check_iR(vname, rval, rmin, rmax)

    character(*), intent(in) :: vname ! variable name (possibly with arbitrary comments)
    real(R_),     intent(in) :: rval  ! value to be checked
    real(R_),     intent(in), optional :: rmin ! min limit
    real(R_),     intent(in), optional :: rmax ! max limit
    ! - Note: Should be rmin <= rval <= rmax, in a normal case.
    character(256) :: msg

    if (present(rmin)) then
       if (rval < rmin) then
          write (msg, *) trim(vname)//' = ', rval, ': Too small. Should be >= ', rmin
          call err_issue(-1, msg)
       end if
    end if

    if (present(rmax)) then
       if (rval > rmax) then
          write (msg, *) trim(vname)//' = ', rval, ': Too large. Should be <= ', rmax
          call err_issue(1, msg)
       end if
    end if

  end subroutine check_iR


  !+
  ! Check whether a real scalar variable has a valid numeric value
  !-
  subroutine check_valid(vname, rval)

    character(*), intent(in) :: vname ! variable name (possibly with arbitrary comments)
    real(R_),     intent(in) :: rval  ! value to be checked
    character(256) :: msg

    if (rval >= -RLRG_ .and. rval <= RLRG_) then
    else
       write (msg, *) trim(vname)//' = ', rval, ': Seems to be a invalid numeric value.'
       call err_issue(-1, msg)
    end if

  end subroutine check_valid


  !+
  ! Convert 4 or less characters to integer(4), with the Little Endian rule
  !-
  function convLtl_AN_to_I(aa) result(j)

    character(*), intent(in) :: aa ! first 4 elements are effective
    integer  :: j  ! result in [-2**31, 2**31-1]
    integer, parameter :: I8 = 256
    integer, parameter :: I16 = I8 * I8
    integer, parameter :: I24 = I8 * I8 * I8
    integer  :: n, i, jj(4)

    ! 4 or less integer values in [-128, 127]
    n = min(4, len(aa))
    do i = 1, n
       jj(i) = iachar(aa(i:i))
    end do

    ! Conversion
    if (n == 4) then
       if (jj(4) <= 127) then
          j = jj(4) * I24
       else
          j = IMIN_ + (jj(4) - 128) * I24
       end if
    else
       j = 0
    end if
    if (n >= 3) j = j + jj(3) * I16
    if (n >= 2) j = j + jj(2) * I8
    j = j + jj(1)

  end function convLtl_AN_to_I


  !+
  ! Convert a character string to an integer vector, with the Little Endian rule
  !-
  function convLtl_AN_to_1I(str) result(ivec)

    character(*), intent(in) :: str     ! source
    integer :: ivec((len(str) - 1) / 4 + 1) ! destination
    integer :: i, n, ics, ice

    n = (len(str) - 1) / 4 + 1
    do i = 1, n
       ics = 4 * i - 3
       ice = min(len(str), 4 * i)
       ivec(i) = convLtl_AN_to_I(str(ics:ice))
    end do

  end function convLtl_AN_to_1I


  !+
  ! Add integer scalars in a scalar, cyclically in the range [-2**31, 2**31-1]
  !-
  function cyclicAdd_I(jdat1, jdat2) result(jdat3)

    integer,  intent(in) :: jdat1     ! integer scalar 1
    integer,  intent(in) :: jdat2     ! integer scalar 2
    integer  :: jdat3     ! output scalar

    if (jdat2 >= 0) then
       if (jdat1 <= IMAX_ - jdat2) then
          jdat3 = jdat1 + jdat2
       else
          jdat3 = (jdat1 - IMAX_ - 1) + (jdat2 - IMAX_ - 1)
       end if
    else
       if (jdat1 >= -IMAX_ - jdat2 - 1) then
          jdat3 = jdat1 + jdat2
       else
          jdat3 = (jdat1 + IMAX_ + 1) + (jdat2 + IMAX_ + 1)
       end if
    end if

  end function cyclicAdd_I


  !+
  ! Set unit index for standard error output
  !-
  subroutine err_setUnit(iue) 

    integer, intent(in) :: iue
    Err_iue = iue

  end subroutine err_setUnit


  !+
  ! Set an action to do when an error is detected
  !-
  subroutine err_setAct(mact) 

    integer, intent(in) :: mact
    Err_mact = mact

  end subroutine err_setAct


  !+
  ! Check status and possibly issue an error message
  !-
  subroutine err_issue(is, msg)

    integer, intent(in) :: is  ! status code (0 = normal, others = errors)
    character(len=*), intent(in), optional :: msg   ! error message

    if (is /= 0 .and. Err_mact >= 0) then
       if (present(msg)) call writeStr(Err_iue, msg)
       write (Err_iue, *) 'Error code = ', is
       if (Err_mact == 1) stop
    end if

  end subroutine err_issue


  !+
  ! Check status after opening a file
  !-
  subroutine err_open(is, iu, msg)

    integer,      intent(in) :: is  ! status code (0 = normal, others = errors)
    integer,      intent(in) :: iu  ! unit index of the file tried to be opened
    character(*), intent(in), optional :: msg   ! error message
    logical :: op
    character(256) :: fname

    if (is /= 0 .and. Err_mact >= 0) then
       if (present(msg)) call writeStr(Err_iue, msg)
       write (Err_iue, *) 'Error code = ', is, ' : Failed to open.'
       inquire (unit=iu, opened=op, name=fname)
       if (op) then
          write (Err_iue, *) 'File index = ', iu, &
               & '. Already opened, being linked to ', trim(fname)
       else
          write (Err_iue, *) 'File index = ', iu
       end if
       if (Err_mact == 1) stop
    end if

  end subroutine err_open


  !+
  ! Check status after reading from a file
  !-
  subroutine err_read(is, iu, msg)

    integer,      intent(in) :: is  ! status code (0 = normal, others = errors)
    integer,      intent(in) :: iu  ! unit index of the file
    character(*), intent(in), optional :: msg   ! error message
    logical :: op
    character(256) :: fname

    if (is /= 0 .and. Err_mact >= 0) then
       if (present(msg)) call writeStr(Err_iue, msg)
       if (is < 0) then
          write (Err_iue, *) 'Error, code = ', is, ' : Failed to read in. EOF was found.'
       else
          write (Err_iue, *) 'Error, code = ', is, ' : Failed to read in.'
       end if
       inquire (unit=iu, opened=op, name=fname)
       if (op) then
          write (Err_iue, *) 'File = ', trim(fname)
       else
          write (Err_iue, *) 'File index = ', iu, '. Not yet opened.'
       end if
       if (Err_mact == 1) stop
    end if

  end subroutine err_read


  !+
  ! Check status after writing to a file
  !-
  subroutine err_write(is, iu, msg)

    integer,      intent(in) :: is  ! status code (0 = normal, others = errors)
    integer,      intent(in) :: iu  ! unit index of the file
    character(*), intent(in), optional :: msg   ! error message
    logical :: op
    character(256) :: fname

    if (is /= 0 .and. Err_mact >= 0) then
       if (present(msg)) call writeStr(Err_iue, msg)
       if (is < 0) then
          write (Err_iue, *) 'Error, code = ', is, ' : Failed to writing out. EOF was found.'
       else
          write (Err_iue, *) 'Error, code = ', is, ' : Failed to writing out.'
       end if
       inquire (unit=iu, opened=op, name=fname)
       if (op) then
          write (Err_iue, *) 'File = ', trim(fname)
       else
          write (Err_iue, *) 'File index = ', iu, '. Not yet opened.'
       end if
       if (Err_mact == 1) stop
    end if

  end subroutine err_write


  !+
  ! Returns a file name without path
  !-
  function fileName_AN(fullname, sep) result(fname)

    character(*), intent(in) :: fullname ! full file name possibly with a path
    character(1), intent(in), optional :: sep ! directory name separator (usually '/')
    character(len(fullname)) :: fname ! result, file name without path
    character(1) :: sep1
    integer   :: i, nf, np

    sep1 = '/'
    if (present(sep)) sep1 = sep
    nf = len_trim(fullname)
    np = 0
    do i = nf, 1, -1
       if (fullname(i:i) == sep1) then ! find the last separator
          np = i
          exit
       end if
    end do
    if (np >= 1) then ! result
       fname = fullname(np+1:)
    else
       fname = fullname
    end if

  end function fileName_AN


  !+
  ! Returns a file name path
  !-
  function filePath_AN(fname, sep) result(pname)

    character(*), intent(in) :: fname ! file name
    character(1), intent(in), optional :: sep ! directory name separator (usually '/')
    character(len(fname)) :: pname ! result, file path name
    character(1) :: sep1
    integer   :: i, nf, np

    sep1 = '/'
    if (present(sep)) sep1 = sep
    nf = len_trim(fname)
    np = 0
    do i = nf, 1, -1
       if (fname(i:i) == sep1) then ! find the last separator
          np = i
          exit
       end if
    end do
    if (np >= 1) then ! result
       pname = fname(1:np)
    else ! no path name
       pname = '.'
    end if

  end function filePath_AN


  !+
  ! Unit index freely available for a file
  !-
  function freeUnit_I(iu0) result(iu)

    integer, intent(in), optional :: iu0 ! a unit index as a candidate
    integer :: iu ! a unit index that can be used
    integer :: iu00
    logical :: op

    iu00 = 10
    if (present(iu0)) iu00 = max(10, iu0)
    op = .true.
    do iu = iu00, 99, 1
       inquire (unit=iu, opened=op)
       if (.not. op) exit
    end do
    if (op) then
       do iu = 10, iu00 - 1
          inquire (unit=iu, opened=op)
          if (.not. op) exit
       end do
       if (op) call err_issue(1, 'freeUnit_I: No unit index is found.')
    end if

  end function freeUnit_I


  !+
  ! Convert an integer variable to a character string
  !-
  function num2str_AN(i, fmt) result(str)

    integer, intent(in) :: i
    character(*), intent(in), optional :: fmt
    character(12) :: str
    character(24) :: tmpstr

    tmpstr = ' '
    if (present(fmt)) then
       write (tmpstr, fmt) i
    else
       write (tmpstr, *) i
    end if
    str = adjustl(tmpstr)

  end function num2str_AN


  !+
  ! A non-zero value with an absolute value larger than threshold
  !-
  function nonZero_R(val, sml) result(res)

    real(R_), intent(in) :: val ! input value
    real(R_), intent(in) :: sml ! threshold
    real(R_) :: res ! output value

    if (val < -sml .or. val > sml) then
       res = val
    else if (val >= 0.0_R_) then
       res = sml
    else
       res = -sml
    end if

  end function nonZero_R


  !+
  ! Return # of digits of integer variable
  !-
  function numDigits_I(i) result(n)

    integer, intent(in) :: i
    integer :: n
    integer :: ii, id, im

    ii = abs(i)
    im = 10
    do id = 1, 19
       if (ii < im) exit
       im = im * 10
    end do
    n = id

  end function numDigits_I


  !+
  ! Read unknown length character string and count its length
  !-
  subroutine read1word(iu, str)

    integer, intent(in) :: iu
    character(*), intent(out) :: str
    integer :: ios, ln

    read (iu, '(a)', iostat=ios) str
    call err_read(ios, iu)
    str = adjustl(str)
    ln = index(str, ' ')
    str = str(1:ln)

  end subroutine read1word


  !+
  ! Shuffle an integer vector, using an index table for rearranging
  !-
  subroutine shuffle_m1I(jdat, jidx, is, ie, iw)

    integer,  intent(inout) :: jdat(:)     ! integer vector
    integer,  intent(in)    :: jidx(:)     ! index table for rearranging 
    integer,  intent(in)    :: is, ie, iw  ! start, end, & width of grid index for jidx
    integer  :: i1, i2, ii, jtmp

    ii = 0
    do i1 = is, ie, iw
       ii = ii + 1
       i2 = jidx(ii) ! modified, Dec 2009
       jtmp     = jdat(i1)
       jdat(i1) = jdat(i2)
       jdat(i2) = jtmp
    end do

  end subroutine shuffle_m1I


  !+
  ! Shuffle an integer vector, using an index table for rearranging
  !-
  subroutine shuffle_m1A(adat, jidx, is, ie, iw)

    character(1), intent(inout) :: adat(:) ! character data vector
    integer,  intent(in) :: jidx(:)        ! index table for rearranging 
    integer,  intent(in) :: is, ie, iw     ! start, end, & width of grid index for jidx
    character(1) :: atmp
    integer  :: i1, i2, ii

    if (iw < 0) then
       ii = -(is - ie - iw) / iw + 1
    else
       ii = 0
    end if
    do i1 = is, ie, iw
       if (iw < 0) then
          ii = ii - 1
       else
          ii = ii + 1
       end if
       i2 = jidx(ii) ! modified, Dec 2009
       atmp     = adat(i1)
       adat(i1) = adat(i2)
       adat(i2) = atmp
    end do

  end subroutine shuffle_m1A


  !+
  ! Sort a vector by the heap method
  !-
  subroutine sort_heap(dat, nn)

    real(R_), intent(inout) :: dat(:)
    integer,  intent(in), optional :: nn
    integer  :: n, i, ichild, iparent, j
    real(R_) :: tmp

    n = size(dat)
    if (present(nn)) n = nn

    do i = 2, n
       tmp = dat(i)
       j = i
       do
          if (j <= 1) exit
          iparent = j / 2
          if (dat(iparent) >= tmp) exit
          dat(j) = dat(iparent)
          j = iparent
       end do
       dat(j) = tmp
    end do

    do i = n - 1, 1, -1
       call swap_mR(dat(1), dat(i + 1))
       j = 1
       ichild = 2
       tmp = dat(j)
       do
          if (ichild > i) exit
          if (ichild < i .and. dat(ichild + 1) > dat(ichild)) ichild = ichild + 1
          if (tmp >= dat(ichild)) exit
          dat(j) = dat(ichild)
          j = ichild
          ichild = j * 2
       end do
       dat(j) = tmp
    end do

  end subroutine sort_heap


  !+
  ! Sort a vector by the selection method
  !-
  subroutine sort_selec(dat, nn)

    real(R_), intent(inout) :: dat(:)
    integer,  intent(in), optional :: nn
    integer :: n, i0, i1, j

    n = size(dat)
    if (present(nn)) n = nn

    do i0 = 1, n - 1
       j = i0
       do i1 = i0 + 1, n
          if (dat(i1) < dat(j)) j = i1 
       end do
       call swap_mR(dat(i0), dat(j))
    end do

  end subroutine sort_selec


  !+
  ! Convert a character string to an integer variable
  !-
  function str2num_I(str) result(i)

    character(*), intent(in) :: str
    integer :: i
    integer :: ios

    read (str, *, iostat=ios) i
    call err_issue(ios, 'str2num_I: Could not find an integer value from the argument, '//trim(str))

  end function str2num_I


  !+
  ! Convert a character string to an integer value, by skipping the first word
  !// str2num_skip_I('  variableA1 51 99 variable A1') will be = 51
  !-
  function str2num_skip_I(str) result(i)

    character(*), intent(in) :: str
    integer :: i
    integer :: ns, ios

    ns = len(str) - len(adjustl(str)) ! (location of the first word) - 1
    ns = index(str(ns+1:), ' ') + ns  ! location of the first word end
    read (str(ns+1:), *, iostat=ios) i
    call err_issue(ios, 'str2num_skip_I: No integer value is found.')

  end function str2num_skip_I


  !+
  ! Extract a character string that appears after some character string
  !-
  function str_afterStr_AN(str1, str2) result(newstr)

    character(*), intent(in) :: str1
    character(*), intent(in) :: str2
    character(len(str1) - len(str2)) :: newstr
    integer :: ip

    ip = index(str1, str2)
    if (ip /= 0) then
       ip = ip + len(str2)
       if (ip <= len(str1)) then
          newstr = str1(ip:)
       else
          newstr = ' '
       end if
    else
       newstr = ' '
    end if

  end function str_afterStr_AN


  !+
  ! Shift a character string, extracting a character string that contains the second and subsequent words
  !// str_shift_AN('  variableA1 51 99 variable A1') will be = '51 99 variable A1'
  !-
  function str_shift_AN(str, n) result(res)

    character(*), intent(in) :: str ! source string
    integer,      intent(in) :: n   ! # of shift operations
    character(len(str)) :: res
    integer :: ns, i
    
    res = str
    do i = 1, n
       ns = len_trim(res) - len_trim(adjustl(res)) ! (location of the first word) - 1
       ns = index(res(ns+1:), ' ') + ns  ! location of the first word end
       res = adjustl(res(ns+1:)) ! contains the second and subsequent words
    end do

  end function str_shift_AN


  !+
  ! Split a string into words and get start and end indexes for the word location
  !-
  subroutine str_split(str, nw, spc, ics, ice)

    character(*), intent(in)  :: str     ! source string
    integer,      intent(out) :: nw      ! # of words found
    character(1), optional, intent(in)  :: spc    ! separater character (default = ' ')
    integer,      optional, intent(out) :: ics(:) ! start index for location of each word
    integer,      optional, intent(out) :: ice(:) ! end   index for location of each word
    character(1) :: s
    integer :: id(len(str)+1), ic, nwmax

    ! Initialize
    s = ' '
    if (present(spc)) s = spc
    id(:) = 1 ! not a separater
    do ic = 1, len(str)
       if (str(ic:ic) == s) id(ic) = 0 ! a separater
    end do
    nwmax = len(str)
    if (present(ics)) nwmax = size(ics)
    !print *, id(1:len(str))

    ! Scan the string
    nw = 0
    if (id(1) == 1) then
       nw = nw + 1
       if (nw > nwmax) return
       if (present(ics)) ics(nw) = 1
    end if
    do ic = 2, len(str)+1
       if (id(ic)-id(ic-1) == 1) then ! start of a new word
          !print *, nw, ic
          nw = nw + 1
          if (nw > nwmax) return
          if (present(ics)) ics(nw) = ic
       end if
       if (id(ic)-id(ic-1) == -1 .and. present(ics)) ice(nw) = ic-1 ! end of current word
    end do

  end subroutine str_split


  !+
  ! Swap two scalars
  !-
  subroutine swap_mI(i1, i2)

    integer, intent(inout) :: i1, i2
    integer :: it
    it = i1
    i1 = i2
    i2 = it

  end subroutine swap_mI


  !+
  ! Swap two scalars
  !-
  subroutine swap_mR(v1, v2) 

    real(R_), intent(inout) :: v1, v2
    real(R_) :: tt
    tt = v1
    v1 = v2
    v2 = tt

  end subroutine swap_mR


  !+
  ! Swap two vectors
  !-
  subroutine swap_m1R(v1, v2) 

    real(R_), intent(inout) :: v1(:), v2(:)
    real(R_) :: tt(size(v1))
    tt(:) = v1(:)
    v1(:) = v2(:)
    v2(:) = tt(:)

  end subroutine swap_m1R


  !+
  ! Swap two matrices
  !-
  subroutine swap_m2R(w1, w2) 

    real(R_), intent(inout) :: w1(:,:), w2(:,:)
    real(R_) :: tt(size(w1,1), size(w1,2))
    tt(:,:) = w1(:,:)
    w1(:,:) = w2(:,:)
    w2(:,:) = tt(:,:)

  end subroutine swap_m2R


  !+
  ! Output a character string possibly with formatting
  !-
  subroutine writeStr(iu, str) 

    integer,      intent(in) :: iu  ! output device index
    character(*), intent(in) :: str ! character string
    integer :: ie, is, ln

    ln = len_trim(str)
    if (ln <= 2) then
       write (iu, '(a)') str(1:ln)
    else
       is = 1
       ie = 3
       do ie = 3, ln - 1
          if (str(ie-1:ie) == '\n') then
             if (is <= ie-2) then
                write (iu, '(a)') str(is:ie-2)
             else
                write (iu, *)
             end if
             is = ie + 1
          end if
       end do
       write (iu, '(a)') str(is:ln)
    end if

  end subroutine writeStr

end module hparx_base
