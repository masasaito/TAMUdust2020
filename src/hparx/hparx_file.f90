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
! Library of file I/O utilities
!-
module hparx_file

  use globals, only : R_, RD_, R4_, R8_, I4_, I8_
  use hparx_base
  implicit none
  private

  ! Public
  public :: bin_read_o1R4
  public :: bin_read_o2R4
  public :: bin_read_o3R4
  public :: bin_write_i1R4
  public :: bin_write_i2R4
  public :: bin_write_i3R4
  public :: bin1B_read_o1I
  public :: bin1B_write_i1I
  public :: close_all
  public :: col_read_o1R
  public :: col_read_o1Rx2
  public :: col_read_o1Rx3
  public :: col_read_o1Rx4
  public :: col_write_i1R
  public :: col_write_i1Rx2
  public :: col_write_i1Rx3
  public :: col_write_i1Rx4
  public :: fileRec_init
  public :: fileRec_final
  public :: fileRec_check
  public :: fileRec_set
  public :: fileRec_qry_I
  public :: gradsCtl_read
  public :: gradsCtl_write
  public :: gradsCtl2_read
  public :: readAnnotLine
  public :: open_dir
  public :: open_seq
  public :: recLen_I
  public :: recordLen_I
  public :: row_read_o2R
  public :: row_read_o3R
  public :: row_write_i2R
  public :: row_write_i3R
  public :: skipLines
  public :: skipToSign
  public :: skipToSigns

  ! Private variables: status table for file units
  integer, save, allocatable :: Rec_idx(:) ! table of record indexes for direct-access files
  integer, save :: Rec_nuni = 1000         ! default allocation size for tables
  !   This module is safe to call any procedure without explicit call of the initialization
  !   (fileRec_init). However, if file units of more than the default number are occasionally used,
  !   then call first the initialization procedure with a proper number before calling other procedures.

contains

  !+
  ! Read in 1-D real(R4_) binary data from a direct-access file
  ! - Features usability and instant byte swapping.
  !-
  subroutine bin_read_o1R4(iu, dat, irec, nx, mswap, iostat)

    integer,   intent(in)  :: iu     ! file unit index
    real(R4_), intent(out) :: dat(:) ! data
    integer,   intent(in), optional :: irec   ! record index
    integer,   intent(in), optional :: nx     ! dimension size
    integer,   intent(in), optional :: mswap  ! flag for byte swapping (0=no, 2/4/8=yes)
    integer,   intent(out), optional :: iostat ! I/O status (0:normal, <0:EOF, >0:error)
    !// Default: irec=(auto), nx=size(dat,1), mswap=0
    !// Without iostat option, this code will stop when a read error is detected.
    integer :: irec1, nx1, ios, ix

    ! Initialize
    if (present(irec)) then
       irec1 = irec
    else
       irec1 = fileRec_qry_I(iu) + 1
    end if
    if (present(nx)) then
       nx1 = nx
    else
       nx1 = size(dat, 1)
    end if

    ! Read in
    read (iu, rec=irec1, iostat=ios) (dat(ix), ix = 1, nx1)
    if (present(iostat)) then
       iostat = ios
    else if (ios /= 0) then
       call err_read(ios, iu, 'bin_read_o1R4 : irec='//num2str_AN(irec1))
    end if
    call fileRec_set(iu, irec1) ! save the record index

    ! Byte swap
    if (present(mswap)) then
       if (mswap >= 2) call bswap_m1R4(dat(1:nx1), nx1, mswap)
    end if

  end subroutine bin_read_o1R4


  !+
  ! Read in 2-D real(R4_) binary data from a direct-access file
  ! - Features usability and instant byte swapping.
  !-
  subroutine bin_read_o2R4(iu, dat, mfmt, irec, nx, ny, mswap, iostat)

    integer,   intent(in)  :: iu       ! file unit index
    real(R4_), intent(out) :: dat(:,:) ! data
    integer,   intent(in), optional :: mfmt   ! flag for single record data format (1=1-D, 2=2-D)
    integer,   intent(in), optional :: irec   ! record index
    integer,   intent(in), optional :: nx, ny ! dimension sizes
    integer,   intent(in), optional :: mswap  ! flag for byte swapping (0=no, 2/4/8=yes)
    integer,   intent(out), optional :: iostat ! I/O status (0:normal, <0:EOF, >0:error)
    !// Default: mfmt=2, irec=(auto), nx=size(dat,1), ny=size(dat,2), mswap=0
    !// Without iostat option, this code will stop when a read error is detected.
    integer :: irec1, nx1, ny1, ios, ix, iy, mfmt1

    ! Initialize
    if (present(irec)) then
       irec1 = irec
    else
       irec1 = fileRec_qry_I(iu) + 1
    end if
    if (present(nx)) then
       nx1 = nx
    else
       nx1 = size(dat, 1)
    end if
    if (present(ny)) then
       ny1 = ny
    else
       ny1 = size(dat, 2)
    end if
    mfmt1 = 2
    if (present(mfmt)) mfmt1 = mfmt
    if (present(iostat)) iostat = 0

    ! Read in
    if (mfmt1 == 1) then ! single record = 1-D
       do iy = 1, ny1
          read (iu, rec=irec1, iostat=ios) (dat(ix, iy), ix = 1, nx1)
          if (ios /= 0) then ! on error
             if (present(iostat)) then
                iostat = ios
                exit
             else
                call err_read(ios, iu, 'bin_read_o2R4 : irec='//num2str_AN(irec1))
             end if
          end if
          irec1 = irec1 + 1
       end do
       irec1 = irec1 - 1
    else                ! single record = 2-D
       read (iu, rec=irec1, iostat=ios) ((dat(ix, iy), ix = 1, nx1), iy = 1, ny1)
       if (present(iostat)) then
          iostat = ios
       else if (ios /= 0) then
          call err_read(ios, iu, 'bin_read_o2R4 : irec='//num2str_AN(irec1))
       end if
    end if
    call fileRec_set(iu, irec1) ! save the record index

    ! Byte swap
    if (present(mswap)) then
       if (mswap >= 2) then
          do iy = 1, ny1
             call bswap_m1R4(dat(1:nx1, iy), nx1, mswap)
          end do
       end if
    end if

  end subroutine bin_read_o2R4


  !+
  ! Read in 3-D real(R4_) binary data from a direct-access file
  ! - Features usability and instant byte swapping.
  !-
  subroutine bin_read_o3R4(iu, dat, mfmt, irec, nx, ny, nz, mswap, iostat)

    integer,   intent(in)  :: iu         ! file unit index
    real(R4_), intent(out) :: dat(:,:,:) ! data
    integer,   intent(in), optional :: mfmt ! flag for single record data format (1=1-D, 2=2-D, 3=3-D)
    integer,   intent(in), optional :: irec ! record index
    integer,   intent(in), optional :: nx, ny, nz ! dimension sizes
    integer,   intent(in), optional :: mswap      ! flag for byte swapping (0=no, 2/4/8=yes)
    integer,   intent(out), optional :: iostat ! I/O status (0:normal, <0:EOF, >0:error)
    !// Default: mfmt=3, irec=(auto), nx=size(dat,1), ny=size(dat,2), ny=size(dat,3), mswap=0
    !// Without iostat option, this code will stop when a read error is detected.
    integer :: irec1, nx1, ny1, nz1, ios, ix, iy, iz, mfmt1

    ! Initialize
    if (present(irec)) then
       irec1 = irec
    else
       irec1 = fileRec_qry_I(iu) + 1
    end if
    if (present(nx)) then
       nx1 = nx
    else
       nx1 = size(dat, 1)
    end if
    if (present(ny)) then
       ny1 = ny
    else
       ny1 = size(dat, 2)
    end if
    if (present(nz)) then
       nz1 = nz
    else
       nz1 = size(dat, 3)
    end if
    mfmt1 = 3
    if (present(mfmt)) mfmt1 = mfmt
    if (present(iostat)) iostat = 0

    ! Read in
    if (mfmt1 == 1) then      ! single record = 1-D
       do iz = 1, nz1
          do iy = 1, ny1
             read (iu, rec=irec1, iostat=ios) (dat(ix, iy, iz), ix = 1, nx1)
             if (ios /= 0) then ! on error
                if (present(iostat)) then
                   iostat = ios
                   exit
                else
                   call err_read(ios, iu, 'bin_read_o3R4 : irec='//num2str_AN(irec1))
                end if
             end if
             irec1 = irec1 + 1
          end do
       end do
       irec1 = irec1 - 1

    else if (mfmt1 == 2) then ! single record = 2-D
       do iz = 1, nz1
          read (iu, rec=irec1, iostat=ios) ((dat(ix, iy, iz), ix = 1, nx1), iy = 1, ny1)
          if (ios /= 0) then ! on error
             if (present(iostat)) then
                iostat = ios
                exit
             else
                call err_read(ios, iu, 'bin_read_o3R4 : irec='//num2str_AN(irec1))
             end if
          end if
          irec1 = irec1 + 1
       end do
       irec1 = irec1 - 1

    else                     ! single record = 3-D
       read (iu, rec=irec1, iostat=ios) (((dat(ix, iy, iz), ix = 1, nx1), iy = 1, ny1), iz = 1, nz1)
       if (present(iostat)) then
          iostat = ios
       else if (ios /= 0) then
          call err_read(ios, iu, 'bin_read_o3R4 : irec='//num2str_AN(irec1))
       end if
    end if
    call fileRec_set(iu, irec1) ! save the record index

    ! Byte swap
    if (present(mswap)) then
       if (mswap >= 2) then
          do iz = 1, nz1
             do iy = 1, ny1
                call bswap_m1R4(dat(1:nx1, iy, iz), nx1, mswap)
             end do
          end do
       end if
    end if

  end subroutine bin_read_o3R4


  !+
  ! Write out 1-D real(R4_) binary data to a direct-access file
  ! - Features usability.
  !-
  subroutine bin_write_i1R4(iu, dat, irec, nx)

    integer,   intent(in) :: iu     ! file unit index
    real(R4_), intent(in) :: dat(:) ! data
    integer,   intent(in), optional :: irec   ! record index
    integer,   intent(in), optional :: nx     ! dimension size
    integer :: irec1, nx1, ios, ix

    ! Initialize
    if (present(irec)) then
       irec1 = irec
    else
       irec1 = fileRec_qry_I(iu) + 1
    end if
    if (present(nx)) then
       nx1 = nx
    else
       nx1 = size(dat, 1)
    end if

    ! Write out
    write (iu, rec=irec1, iostat=ios) (dat(ix), ix = 1, nx1)
    if (ios /= 0) call err_write(ios, iu, 'bin_write_i1R4 : irec='//num2str_AN(irec1))
    call fileRec_set(iu, irec1) ! save the record index

  end subroutine bin_write_i1R4


  !+
  ! Write out 2-D real(R4_) binary data to a direct-access file
  ! - Features usability.
  !-
  subroutine bin_write_i2R4(iu, dat, mfmt, irec, nx, ny)

    integer,   intent(in) :: iu         ! file unit index
    real(R4_), intent(in) :: dat(:,:)   ! data
    integer,   intent(in), optional :: mfmt ! flag for single record data format (1=1-D, 2=2-D)
    integer,   intent(in), optional :: irec ! record index
    integer,   intent(in), optional :: nx, ny   ! dimension sizes
    !// Default: mfmt=2, irec=(auto), nx=size(dat,1), ny=size(dat,2), mswap=0
    integer :: irec1, nx1, ny1, ios, ix, iy, mfmt1

    ! Initialize
    if (present(irec)) then
       irec1 = irec
    else
       irec1 = fileRec_qry_I(iu) + 1
    end if
    if (present(nx)) then
       nx1 = nx
    else
       nx1 = size(dat, 1)
    end if
    if (present(ny)) then
       ny1 = ny
    else
       ny1 = size(dat, 2)
    end if
    mfmt1 = 2
    if (present(mfmt)) mfmt1 = mfmt

    ! Write out
    if (mfmt1 == 1) then    ! single record = 1-D
       do iy = 1, ny1
          write (iu, rec=irec1, iostat=ios) (dat(ix, iy), ix = 1, nx1)
          if (ios /= 0) call err_write(ios, iu, 'bin_write_i3R4 : irec='//num2str_AN(irec1))
          irec1 = irec1 + 1
       end do
       irec1 = irec1 - 1
    else                   ! single record = 2-D
       write (iu, rec=irec1, iostat=ios) ((dat(ix, iy), ix = 1, nx1), iy = 1, ny1)
       if (ios /= 0) call err_write(ios, iu, 'bin_write_i3R4 : irec='//num2str_AN(irec1))
    end if
    call fileRec_set(iu, irec1) ! save the record index

  end subroutine bin_write_i2R4


  !+
  ! Write out 3-D real(R4_) binary data to a direct-access file
  ! - Features usability.
  !-
  subroutine bin_write_i3R4(iu, dat, mfmt, irec, nx, ny, nz)

    integer,   intent(in) :: iu         ! file unit index
    real(R4_), intent(in) :: dat(:,:,:) ! data
    integer,   intent(in), optional :: mfmt ! flag for single record data format (1=1-D, 2=2-D, 3=3-D)
    integer,   intent(in), optional :: irec ! record index
    integer,   intent(in), optional :: nx, ny, nz ! dimension sizes
    !// Default: mfmt=3, irec=(auto), nx=size(dat,1), ny=size(dat,2), nz=size(dat,3), mswap=0
    integer :: irec1, nx1, ny1, nz1, ios, ix, iy, iz, mfmt1

    ! Initialize
    if (present(irec)) then
       irec1 = irec
    else
       irec1 = fileRec_qry_I(iu) + 1
    end if
    if (present(nx)) then
       nx1 = nx
    else
       nx1 = size(dat, 1)
    end if
    if (present(ny)) then
       ny1 = ny
    else
       ny1 = size(dat, 2)
    end if
    if (present(nz)) then
       nz1 = nz
    else
       nz1 = size(dat, 3)
    end if
    mfmt1 = 3
    if (present(mfmt)) mfmt1 = mfmt

    ! Write out
    if (mfmt1 == 1) then      ! single record = 1-D
       do iz = 1, nz1
          do iy = 1, ny1
             write (iu, rec=irec1, iostat=ios) (dat(ix, iy, iz), ix = 1, nx1)
             if (ios /= 0) call err_write(ios, iu, 'bin_write_i3R4 : irec='//num2str_AN(irec1))
             irec1 = irec1 + 1
          end do
       end do
       irec1 = irec1 - 1
    else if (mfmt1 == 2) then ! single record = 2-D
       do iz = 1, nz1
          write (iu, rec=irec1, iostat=ios) ((dat(ix, iy, iz), ix = 1, nx1), iy = 1, ny1)
          if (ios /= 0) call err_write(ios, iu, 'bin_write_i3R4 : irec='//num2str_AN(irec1))
          irec1 = irec1 + 1
       end do
       irec1 = irec1 - 1
    else                     ! single record = 3-D
       write (iu, rec=irec1, iostat=ios) (((dat(ix, iy, iz), ix = 1, nx1), iy = 1, ny1), iz = 1, nz1)
       if (ios /= 0) call err_write(ios, iu, 'bin_write_i3R4 : irec='//num2str_AN(irec1))
    end if
    call fileRec_set(iu, irec1) ! save the record index

  end subroutine bin_write_i3R4


  !+
  ! Read in 1-byte binary data stream and convert to integer vector
  !-
  subroutine bin1B_read_o1I(iu, jdat, irec, nx)

    integer,  intent(in)  :: iu      ! file unit index
    integer,  intent(out) :: jdat(:) ! data (0-255 integer)
    integer,  intent(in), optional :: irec ! record index
    integer,  intent(in), optional :: nx   ! # of vector elements to be read
    character(1) :: abuf(size(jdat)) !(AUTO)
    integer :: irec1, nx1, ios, ix

    ! Initialize
    if (present(irec)) then
       irec1 = irec
    else
       irec1 = fileRec_qry_I(iu) + 1
    end if
    if (present(nx)) then
       nx1 = nx
    else
       nx1 = size(jdat, 1)
    end if

    ! Input
    read (iu, rec=irec1, iostat=ios) (abuf(ix), ix = 1, nx1)
    call err_read(ios, iu, 'bin1B_read_o1I : irec='//num2str_AN(irec1))
    call fileRec_set(iu, irec1) ! save the record index
    jdat(:) = ichar(abuf(:)) ! convert

  end subroutine bin1B_read_o1I


  !+
  ! Write out integer vector data to 1-byte binary data stream
  !-
  subroutine bin1B_write_i1I(iu, jdat, irec, nx)

    integer,  intent(in) :: iu      ! file unit index
    integer,  intent(in) :: jdat(:) ! data (0-255 integer)
    integer,  intent(in), optional :: irec ! record index
    integer,  intent(in), optional :: nx   ! # of vector elements to be written
    character(1) :: abuf(size(jdat)) !(AUTO)
    integer :: irec1, nx1, ios, ix

    ! Initialize
    if (present(irec)) then
       irec1 = irec
    else
       irec1 = fileRec_qry_I(iu) + 1
    end if
    if (present(nx)) then
       nx1 = nx
    else
       nx1 = size(jdat, 1)
    end if

    ! Write out
    abuf(:) = char(jdat(:))
    write (iu, rec=irec1, iostat=ios) (abuf(ix), ix = 1, nx1)
    if (ios /= 0) call err_read(ios, iu, 'bin_write_i1R4 : irec='//num2str_AN(irec1))
    call fileRec_set(iu, irec1) ! save the record index

  end subroutine bin1B_write_i1I


  !+
  ! Close files
  !-
  subroutine close_all(iiu)

    integer, intent(in) :: iiu(:)  ! file unit index
    integer :: i
    do i = 1, size(iiu)
       close(iiu(i))
    end do
    
  end subroutine close_all


  !+
  ! Read in 1 column from file
  !  Features automatic header skipping, free choice of read column, and automatic 
  !  detection of data sequence size 
  !-
  subroutine col_read_o1R(iu, nx, dat1, jc, nxread)

    integer,  intent(in)  :: iu  ! file unit index
    integer,  intent(in)  :: nx  ! # of rows to be read (0=unknown)
    real(R_), intent(out) :: dat1(:)  ! data read
    integer,  intent(in),  optional :: jc  ! indexes of columns to be read
    integer,  intent(out), optional :: nxread ! # of rows actualy read (if nxread < nx, error!)
    real(R_) :: wrk(100)
    integer :: jctab(1), nc, ic, ios, ix, nxqry, nx1

    ! Setup
    if (present(jc)) then
       jctab(1) = jc
    else
       jctab(1) = 1
    end if
    nc = jctab(1)      ! # of buffers to be read
    if (nx >= 1) then
       nxqry = nx      ! # of query data rows
    else
       nxqry = size(dat1) ! unknown # of data
    end if
    nx1 = 0

    ! The first row
    do
       read (iu, *, iostat=ios) (wrk(ic), ic = 1, nc)
       if (ios <= 0) exit
    end do
    if (ios < 0) then
       call err_read(ios, iu, 'col_read_o1R : No data found.')
       if (present(nxread)) nxread = nx1
       return
    end if
    dat1(1) = wrk(jctab(1))
    nx1 = 1

    ! Read in subsequent rows
    do ix = 2, nxqry
       read (iu, *, iostat=ios) (wrk(ic), ic = 1, nc)
       if (ios /= 0) exit
       dat1(ix) = wrk(jctab(1))
       nx1 = nx1 + 1
    end do
    if (nx1 < nx) call err_read(ios, iu, 'col_read_o1R : Too few data.')
    if (present(nxread)) nxread = nx1

  end subroutine col_read_o1R


  !+
  ! Read in 2 columns from file
  !  Features automatic header skipping, free choice of read columns, and automatic 
  !  detection of data sequence size 
  !-
  subroutine col_read_o1Rx2(iu, nx, dat1, dat2, jc, nxread)

    integer,  intent(in)  :: iu      ! file unit index
    integer,  intent(in)  :: nx      ! # of rows to be read (0=unknown)
    real(R_), intent(out) :: dat1(:), dat2(:) ! data read
    integer,  intent(in),  optional :: jc(2)  ! indexes of columns to be read
    integer,  intent(out), optional :: nxread ! # of rows actualy read (if nxread < nx, error!)
    real(R_) :: wrk(100)
    integer :: jctab(2), nc, ic, ios, ix, nxqry, nx1

    ! Setup
    if (present(jc)) then
       jctab(1:2) = jc(1:2)
    else
       jctab(1:2) = (/1, 2/)
    end if
    nc = maxval(jctab) ! # of buffers to be read
    if (nx >= 1) then
       nxqry = nx      ! # of query data rows
    else
       nxqry = size(dat1)
    end if
    nx1 = 0

    ! The first row
    do
       read (iu, *, iostat=ios) (wrk(ic), ic = 1, nc)
       if (ios <= 0) exit
    end do
    if (ios < 0) then
       call err_read(ios, iu, 'col_read_o1Rx2 : No data found.')
       if (present(nxread)) nxread = nx1
       return
    end if
    dat1(1) = wrk(jctab(1))
    dat2(1) = wrk(jctab(2))
    nx1 = 1

    ! Read in subsequent rows
    do ix = 2, nxqry
       read (iu, *, iostat=ios) (wrk(ic), ic = 1, nc)
       if (ios /= 0) exit
       dat1(ix) = wrk(jctab(1))
       dat2(ix) = wrk(jctab(2))
       nx1 = nx1 + 1
    end do
    if (nx1 < nx) call err_read(ios, iu, 'col_read_o1Rx2 : Too few data.')
    if (present(nxread)) nxread = nx1

  end subroutine col_read_o1Rx2


  !+
  ! Read in 3 columns from file
  !  Features automatic header skipping, free choice of read columns, and automatic 
  !  detection of data sequence size 
  !-
  subroutine col_read_o1Rx3(iu, nx, dat1, dat2, dat3, jc, nxread)

    integer,  intent(in)  :: iu      ! file unit index
    integer,  intent(in)  :: nx      ! # of rows to be read (0=unknown)
    real(R_), intent(out) :: dat1(:), dat2(:), dat3(:) ! data read
    integer,  intent(in),  optional :: jc(3) ! indexes of columns to be read
    integer,  intent(out), optional :: nxread ! # of rows actualy read (if nxread < nx, error!)
    real(R_) :: wrk(100)
    integer :: jctab(3), nc, ic, ios, ix, nxqry, nx1

    ! Setup
    if (present(jc)) then
       jctab(1:3) = jc(1:3)
    else
       jctab(1:3) = (/1, 2, 3/)
    end if
    nc = maxval(jctab) ! # of buffers to be read
    if (nx >= 1) then
       nxqry = nx      ! # of query data rows
    else
       nxqry = size(dat1)
    end if
    nx1 = 0

    ! The first row
    do
       read (iu, *, iostat=ios) (wrk(ic), ic = 1, nc)
       if (ios <= 0) exit
    end do
    if (ios < 0) then
       call err_read(ios, iu, 'col_read_o1Rx3 : No data found.')
       if (present(nxread)) nxread = nx1
       return
    end if
    dat1(1) = wrk(jctab(1))
    dat2(1) = wrk(jctab(2))
    dat3(1) = wrk(jctab(3))
    nx1 = 1

    ! Read in subsequent rows
    do ix = 2, nxqry
       read (iu, *, iostat=ios) (wrk(ic), ic = 1, nc)
       if (ios /= 0) exit
       dat1(ix) = wrk(jctab(1))
       dat2(ix) = wrk(jctab(2))
       dat3(ix) = wrk(jctab(3))
       nx1 = nx1 + 1
    end do
    if (nx1 < nx) call err_read(ios, iu, 'col_read_o1Rx3 : Too few data.')
    if (present(nxread)) nxread = nx1

  end subroutine col_read_o1Rx3


  !+
  ! Read in 4 columns from file
  !  Features automatic header skipping, free choice of read columns, and automatic 
  !  detection of data sequence size 
  !-
  subroutine col_read_o1Rx4(iu, nx, dat1, dat2, dat3, dat4, jc, nxread)

    integer,  intent(in)  :: iu            ! file unit index
    integer,  intent(in)  :: nx            ! # of rows to be read (0=unknown)
    real(R_), intent(out) :: dat1(:), dat2(:), dat3(:), dat4(:) ! data read
    integer,  intent(in),  optional :: jc(4)  ! indexes of columns to be read
    integer,  intent(out), optional :: nxread ! # of rows actualy read (if nxread < nx, error!)
    real(R_) :: wrk(100)
    integer :: jctab(4), nc, ic, ios, ix, nxqry, nx1

    ! Setup
    if (present(jc)) then
       jctab(1:4) = jc(1:4)
    else
       jctab(1:4) = (/1, 2, 3, 4/)
    end if
    nc = maxval(jctab) ! # of buffers to be read
    if (nx >= 1) then
       nxqry = nx      ! # of query data rows
    else
       nxqry = size(dat1)
    end if
    nx1 = 0

    ! The first row
    do
       read (iu, *, iostat=ios) (wrk(ic), ic = 1, nc)
       if (ios <= 0) exit
    end do
    if (ios < 0) then
       call err_read(ios, iu, 'col_read_o1Rx4 : No data found.')
       return
    end if
    dat1(1) = wrk(jctab(1))
    dat2(1) = wrk(jctab(2))
    dat3(1) = wrk(jctab(3))
    dat4(1) = wrk(jctab(4))
    nx1 = 1

    ! Read in subsequent rows
    do ix = 2, nxqry
       read (iu, *, iostat=ios) (wrk(ic), ic = 1, nc)
       if (ios /= 0) exit
       dat1(ix) = wrk(jctab(1))
       dat2(ix) = wrk(jctab(2))
       dat3(ix) = wrk(jctab(3))
       dat4(ix) = wrk(jctab(4))
       nx1 = nx1 + 1
    end do
    if (nx1 < nx) call err_read(ios, iu, 'col_read_o1Rx4 : Too few data.')
    if (present(nxread)) nxread = nx1

  end subroutine col_read_o1Rx4


  !+
  ! Write out 1 data column to a file
  !  Features optional indexing, automatic detection of data sequence size, 
  !   and user-specified formatting
  !-
  subroutine col_write_i1R(iu, dat1, nx, midx, fmt)

    integer,  intent(in) :: iu        ! file unit index
    real(R_), intent(in) :: dat1(:)   ! data to be written
    integer,      intent(in), optional :: nx   ! # of rows to be written
    integer,      intent(in), optional :: midx ! method for formatting
    character(*), intent(in), optional :: fmt  ! format of writing a single datum value (e.g., 'es12.4')
    character(80) :: fmtstr
    integer :: ios, ix, nx1, midx1

    ! Setup
    nx1 = size(dat1)
    midx1 = 0
    fmtstr = ' '
    if (present(nx)) nx1 = nx
    if (present(midx)) midx1 = midx
    if (present(fmt)) then
       if (midx1 >= 1) fmtstr = '(1i8, 1'//trim(fmt)//')'
       if (midx1 <= 0) fmtstr = '(     1'//trim(fmt)//')'
    end if

    ! Write
    do ix = 1, nx1
       if (midx1 <= 0) then
          if (fmtstr == ' ') write (iu, *,      iostat=ios) dat1(ix)
          if (fmtstr /= ' ') write (iu, fmtstr, iostat=ios) dat1(ix)
       else
          if (fmtstr == ' ') write (iu, *,      iostat=ios) ix, dat1(ix)
          if (fmtstr /= ' ') write (iu, fmtstr, iostat=ios) ix, dat1(ix)
       end if
       if (ios /= 0) exit
    end do
    if (ios /= 0) call err_write(ios, iu, 'col_write_i1R :')

  end subroutine col_write_i1R


  !+
  ! Write out 2 data columns to a file
  !  Features optional indexing, automatic detection of data sequence size, 
  !   and user-specified formatting
  !-
  subroutine col_write_i1Rx2(iu, dat1, dat2, nx, midx, fmt)

    integer,  intent(in) :: iu            ! file unit index
    real(R_), intent(in) :: dat1(:), dat2(:)   ! data to be written
    integer,      intent(in), optional :: nx   ! # of rows to be written
    integer,      intent(in), optional :: midx ! method for formatting
    character(*), intent(in), optional :: fmt  ! format of writing a single datum value (e.g., 'es12.4')
    character(80) :: fmtstr
    integer :: ios, ix, nx1, midx1

    ! Setup
    nx1 = min(size(dat1), size(dat2))
    midx1 = 0
    fmtstr = ' '
    if (present(nx)) nx1 = nx
    if (present(midx)) midx1 = midx
    if (present(fmt)) then
       if (midx1 >= 1) fmtstr = '(1i8, 2'//trim(fmt)//')'
       if (midx1 <= 0) fmtstr = '(     2'//trim(fmt)//')'
    end if

    ! Write
    do ix = 1, nx1
       if (midx1 <= 0) then
          if (fmtstr == ' ') write (iu, *,      iostat=ios) dat1(ix), dat2(ix)
          if (fmtstr /= ' ') write (iu, fmtstr, iostat=ios) dat1(ix), dat2(ix)
       else
          if (fmtstr == ' ') write (iu, *,      iostat=ios) ix, dat1(ix), dat2(ix)
          if (fmtstr /= ' ') write (iu, fmtstr, iostat=ios) ix, dat1(ix), dat2(ix)
       end if
       if (ios /= 0) exit
    end do
    if (ios /= 0) call err_write(ios, iu, 'col_write_i1Rx2 :')

  end subroutine col_write_i1Rx2


  !+
  ! Write out 3 data columns to a file
  !  Features optional indexing, automatic detection of data sequence size, 
  !   and user-specified formatting
  !-
  subroutine col_write_i1Rx3(iu, dat1, dat2, dat3, nx, midx, fmt)

    integer,  intent(in) :: iu            ! file unit index
    real(R_), intent(in) :: dat1(:), dat2(:), dat3(:) ! data to be written
    integer,      intent(in), optional :: nx   ! # of rows to be written
    integer,      intent(in), optional :: midx ! method for formatting
    character(*), intent(in), optional :: fmt  ! format of writing a single datum value (e.g., 'es12.4')
    character(80) :: fmtstr
    integer :: ios, ix, nx1, midx1

    ! Setup
    nx1 = min(size(dat1), size(dat2), size(dat3))
    midx1 = 0
    fmtstr = ' '
    if (present(nx)) nx1 = nx
    if (present(midx)) midx1 = midx
    if (present(fmt)) then
       if (midx1 >= 1) fmtstr = '(1i8, 3'//trim(fmt)//')'
       if (midx1 <= 0) fmtstr = '(     3'//trim(fmt)//')'
    end if

    ! Write
    do ix = 1, nx1
       if (midx1 <= 0) then
          if (fmtstr == ' ') write (iu, *,      iostat=ios) dat1(ix), dat2(ix), dat3(ix)
          if (fmtstr /= ' ') write (iu, fmtstr, iostat=ios) dat1(ix), dat2(ix), dat3(ix)
       else
          if (fmtstr == ' ') write (iu, *,      iostat=ios) ix, dat1(ix), dat2(ix), dat3(ix)
          if (fmtstr /= ' ') write (iu, fmtstr, iostat=ios) ix, dat1(ix), dat2(ix), dat3(ix)
       end if
       if (ios /= 0) exit
    end do
    if (ios /= 0) call err_write(ios, iu, 'col_write_i1Rx3 :')

  end subroutine col_write_i1Rx3


  !+
  ! Write out 4 data columns to a file
  !  Features optional indexing, automatic detection of data sequence size, 
  !   and user-specified formatting
  !-
  subroutine col_write_i1Rx4(iu, dat1, dat2, dat3, dat4, nx, midx, fmt)

    integer,  intent(in)  :: iu            ! file unit index
    real(R_), intent(in) :: dat1(:), dat2(:), dat3(:), dat4(:) ! data to be written
    integer,      intent(in), optional :: nx   ! # of rows to be written
    integer,      intent(in), optional :: midx ! method for formatting
    character(*), intent(in), optional :: fmt  ! format of writing a single datum value (e.g., 'es12.4')
    character(80) :: fmtstr
    integer :: ios, ix, nx1, midx1

    ! Setup
    nx1 = min(size(dat1), size(dat2), size(dat3), size(dat4))
    midx1 = 0
    fmtstr = ' '
    if (present(nx)) nx1 = nx
    if (present(midx)) midx1 = midx
    if (present(fmt)) then
       if (midx1 >= 1) fmtstr = '(1i8, 4'//trim(fmt)//')'
       if (midx1 <= 0) fmtstr = '(     4'//trim(fmt)//')'
    end if

    ! Write
    do ix = 1, nx1
       if (midx1 <= 0) then
          if (fmtstr == ' ') write (iu, *,      iostat=ios) dat1(ix), dat2(ix), dat3(ix), dat4(ix)
          if (fmtstr /= ' ') write (iu, fmtstr, iostat=ios) dat1(ix), dat2(ix), dat3(ix), dat4(ix)
       else
          if (fmtstr == ' ') write (iu, *,      iostat=ios) ix, dat1(ix), dat2(ix), dat3(ix), dat4(ix)
          if (fmtstr /= ' ') write (iu, fmtstr, iostat=ios) ix, dat1(ix), dat2(ix), dat3(ix), dat4(ix)
       end if
       if (ios /= 0) exit
    end do
    if (ios /= 0) call err_write(ios, iu, 'col_write_i1Rx4 :')

  end subroutine col_write_i1Rx4


  !+
  ! Initialize the record-index table
  !-
  subroutine fileRec_init(nuni) 

    integer, intent(in), optional :: nuni ! # of file units
    if (present(nuni)) Rec_nuni = nuni
    if (allocated(Rec_idx)) call fileRec_final()
    allocate (Rec_idx(Rec_nuni))
    Rec_idx(:) = 0

  end subroutine fileRec_init


  !+
  ! Finalize the record-index table
  !-
  subroutine fileRec_final() 

    deallocate (Rec_idx)

  end subroutine fileRec_final


  !+
  ! Check allocation status for the record-index table
  !-
  subroutine fileRec_check(iuni)

    integer, intent(in) :: iuni  ! unit index of the query

    if (.not. allocated(Rec_idx)) call fileRec_init()
    if (iuni > Rec_nuni) then
       call err_setAct(1) ! to stop
       call err_issue(1, 'fileRec_check : Too large file unit index.' &
            & //'Use more large table size for fileRec_init().')
    end if

  end subroutine fileRec_check


  !+
  ! Set a record index at a specific location in the table
  !-
  subroutine fileRec_set(iuni, irec)

    integer, intent(in) :: iuni  ! unit index of the query
    integer, intent(in), optional :: irec  ! index of record just accessed

    call fileRec_check(iuni)
    if (present(irec)) then
       Rec_idx(iuni) = irec
    else
       Rec_idx(iuni) = 0
    end if

  end subroutine fileRec_set


  !+
  ! Returns a record index of a file unit
  !-
  function fileRec_qry_I(iuni) result(irec) 

    integer, intent(in) :: iuni  ! unit index of the query
    integer :: irec
    call fileRec_check(iuni)
    irec = Rec_idx(iuni)

  end function fileRec_qry_I


  !+
  ! Read GrADS control file
  !-
  subroutine gradsCtl_read(iu, infile, nx, ny, nz, nt, nvar, varstr)

    integer,      intent(in)  :: iu
    character(*), intent(out) :: infile
    integer,      intent(out) :: nx, ny, nz, nt
    integer,      intent(out) :: nvar
    character(*), intent(out) :: varstr(:)
    character(1024) :: bufstr
    integer :: ip, ivar

    call readAnnotLine(iu, 'DSET', bufstr)
    ip = index(bufstr, '^')
    if (ip >= 1) then
       infile = bufstr(ip+1:)
    else
       infile = adjustl(bufstr)
    end if
    call readAnnotLine(iu, 'XDEF', bufstr)
    read (bufstr, *, err=1) nx
    call readAnnotLine(iu, 'YDEF', bufstr)
    read (bufstr, *, err=1) ny
    call readAnnotLine(iu, 'ZDEF', bufstr)
    read (bufstr, *, err=1) nz
    call readAnnotLine(iu, 'TDEF', bufstr)
    read (bufstr, *, err=1) nt
    call readAnnotLine(iu, 'VARS', bufstr)
    read (bufstr, *, err=1) nvar
    do ivar = 1, nvar
       read (iu, '(a)', err=1) varstr(ivar)
    end do

    return
1   call err_read(1, iu, 'gradsCtl_read')

  end subroutine gradsCtl_read


  !+
  ! Read a GrADS control file, enhanced version
  !-
  subroutine gradsCtl2_read(iu, infile, nx, ny, nz, nt, tdef, nvar, var_name, var_nz, var_title)

    integer,      intent(in)  :: iu             ! file unit index
    character(*), intent(out) :: infile         ! file name
    integer,      intent(out) :: nx, ny, nz, nt ! # of X/Y/Z/T grids
    character(*), intent(out) :: tdef           ! T definition (e.g., "linear 02:30Z18DEC2014 10mn")
    integer,      intent(out) :: nvar           ! # of variables
    character(*), intent(out) :: var_name(:)    ! variable names
    integer,      intent(out) :: var_nz(:)      ! NZ for each variable
    character(*), intent(out) :: var_title(:)   ! variable titles
    character(1024) :: bufstr
    integer :: ip, ivar, ics(5), ice(5), nw

    call readAnnotLine(iu, 'DSET', bufstr)
    ip = index(bufstr, '^')
    if (ip >= 1) then
       infile = bufstr(ip+1:)
    else
       infile = adjustl(bufstr)
    end if
    call readAnnotLine(iu, 'XDEF', bufstr)
    read (bufstr, *, err=1, end=1) nx
    call readAnnotLine(iu, 'YDEF', bufstr)
    read (bufstr, *, err=1, end=1) ny
    call readAnnotLine(iu, 'ZDEF', bufstr)
    read (bufstr, *, err=1, end=1) nz
    call readAnnotLine(iu, 'TDEF', bufstr)
    read (bufstr, *, err=1, end=1) nt
    tdef = str_shift_AN(bufstr, 1)
    call readAnnotLine(iu, 'VARS', bufstr)
    !print *, nx, ny, nz, nt
    read (bufstr, *, err=1, end=1) nvar
    do ivar = 1, nvar
       read (iu, '(a)', err=1, end=1) bufstr
       call str_split(bufstr, nw, ics=ics, ice=ice) ! split into words
       var_name(ivar)  = bufstr(ics(1):ice(1))  ! variable name is from the 1st word
       var_title(ivar) = bufstr(ics(3):)        ! variable title is from the 3rd word
       read(bufstr(ics(2):ice(2)), *, err=1, end=1) var_nz(ivar) ! NZ is from the 2nd word
    end do

    return
1   call err_read(1, iu, 'gradsCtl2_read')

  end subroutine gradsCtl2_read


  !+
  ! Write GrADS control file
  !-
  subroutine gradsCtl_write(iu, outfile, nx, ny, nz, nt, nvar, varstr, tdef, undef)

    integer,      intent(in) :: iu             ! file unit index
    character(*), intent(in) :: outfile        ! output data file
    integer,      intent(in) :: nx, ny, nz, nt ! # of X, Y, Z, T grid points
    integer,      intent(in) :: nvar           ! # of variables
    character(*), intent(in) :: varstr(:)      ! variable definition
    character(*), optional, intent(in) :: tdef ! T definition (e.g., "linear 02:30Z18DEC2014 10mn")
    real(R_),     optional, intent(in) :: undef ! value for undefined data
    integer :: ivar

    write (iu, '(a)') 'DSET '//trim(adjustl(outfile))
    write (iu, '(a)') '*OPTIONS LITTLE_ENDIAN'
    write (iu, '(a)') '*OPTIONS BIG_ENDIAN'
    write (iu, '(a)') '*OPTIONS template'
    write (iu, '(a)') 'TITLE '//trim(outfile)
    if (present(undef)) then
       write (iu, '(a, 1es14.6)') 'UNDEF ', undef
    else
       write (iu, '(a)') 'UNDEF 1.0e+37'
    endif
    write (iu, '(a, 1i9, a)') 'XDEF ', nx, ' LINEAR 1 1'
    write (iu, '(a, 1i9, a)') 'YDEF ', ny, ' LINEAR 1 1'
    write (iu, '(a, 1i9, a)') 'ZDEF ', nz, ' LINEAR 1 1'
    if (present(tdef)) then
       write (iu, '(a, 1i9, a)') 'TDEF ', nt, ' '//tdef
    else
       write (iu, '(a, 1i9, a)') 'TDEF ', nt, ' LINEAR 00:00Z00JAN0001 1YR'
    end if
    write (iu, '(a, 1i9)') 'VARS ', nvar
    do ivar = 1, nvar
       write (iu, '(a)') trim(varstr(ivar))
    end do
    write (iu, '(a)') 'ENDVARS'

  end subroutine gradsCtl_write


  !+
  ! Open a file as direct-access, unformatted device
  !-
  subroutine open_dir(iu, fname, nrec, stat)

    integer,      intent(in) :: iu
    character(*), intent(in) :: fname
    integer,      intent(in) :: nrec
    character(*), intent(in), optional :: stat
    integer :: ios

    if (present(stat)) then
       open (iu, file=fname, access='direct', form='unformatted', status=stat, recl=nrec, iostat=ios)
    else
       open (iu, file=fname, access='direct', form='unformatted', recl=nrec, iostat=ios)
    end if
    call err_open(ios, iu, fname)
    call fileRec_set(iu, 0) ! record index is saved in memory
    
  end subroutine open_dir


  !+
  ! Open a file as direct-access, unformatted device
  !-
  subroutine open_seq(iu, fname, stat)

    integer,      intent(in) :: iu
    character(*), intent(in) :: fname
    character(*), intent(in), optional :: stat
    integer :: ios

    if (present(stat)) then
       open (iu, file=fname, status=stat, iostat=ios)
    else
       open (iu, file=fname, iostat=ios)
    end if
    call err_open(ios, iu, fname)
    
  end subroutine open_seq


  !+
  ! Read a line annotated by some words and return a next character string
  !  following the annotation words
  !-
  subroutine readAnnotLine(iu, annot, bufstr)

    integer,      intent(in)  :: iu
    character(*), intent(in)  :: annot
    character(*), intent(out) :: bufstr

    do 
       read (iu, '(a)', err=1, end=2) bufstr
       bufstr = adjustl(bufstr)
       if (index(bufstr, annot) /= 0) exit
    end do
    bufstr = str_afterStr_AN(bufstr, annot)

    return
1   call err_read( 1, iu, 'readAnnotLine')
    return
2   call err_read(-1, iu, 'readAnnotLine')

  end subroutine readAnnotLine


  !+
  ! Record length of a record in direct access, unformatted file: Obsolete version
  !-
  function recLen_I(mtype, mkind, nvar) result(nr)

    integer, intent(in) :: mtype  ! type of variable (1/2/3)
    integer, intent(in) :: mkind  ! kind of variable (0 to 4)
    integer, intent(in) :: nvar   ! # of variables
    integer :: nr ! record length
    integer :: m

    ! Warning
    print *, 'recLen_I: This function is obsolete. Modify the code to use recordLen_I instead.'

    ! Inquire
    if (mtype == 1) then       ! real
       m = 0
       if (mkind == 1) m = 1
       if (mkind == 2) m = 2
       if (mkind == 3) m = 4
       if (mkind == 4) m = 8
       nr = recordLen_I(1, nvar, m)
    else if (mtype == 2) then  ! integer
       m = 0
       if (mkind == 3) m = 4
       if (mkind == 4) m = 8
       nr = recordLen_I(2, nvar, m)
    else                       ! character
       nr = recordLen_I(3, nvar)
    end if

  end function recLen_I


  !+
  ! Record length of a record in direct access, unformatted file
  !
  ! Example: n = recordLen_I(1, 1000, 4) ! for 1000*(1.0_R4_)
  !-
  function recordLen_I(mtype, nvar, mkind) result(nr)

    integer, intent(in) :: mtype ! type of variable (1=real, 2=integer, 3=character)
    integer, intent(in) :: nvar  ! # of data per record to be read/write
    integer, intent(in), optional :: mkind ! kind of variable (1/2/4/8 for real, 4/8 for integer)
    integer :: nr ! record length (system dependent)
    integer :: ivar

    ! Inquire
    if (mtype == 1) then      ! real
       if (present(mkind)) then
          if (mkind == 1) then
             inquire (iolength = nr) (1.0_R_, ivar = 1, nvar)
          else if (mkind == 2) then
             inquire (iolength = nr) (1.0_RD_, ivar = 1, nvar)
          else if (mkind == 4) then
             inquire (iolength = nr) (1.0_R4_, ivar = 1, nvar)
          else if (mkind == 8) then
             inquire (iolength = nr) (1.0_R8_, ivar = 1, nvar)
          else
             inquire (iolength = nr) (1.0, ivar = 1, nvar)
          end if
       else
          inquire (iolength = nr) (1.0, ivar = 1, nvar)
       end if
    else if (mtype == 2) then ! integer
       if (present(mkind)) then
          if (mkind == 4) then
             inquire (iolength = nr) (1_I4_, ivar = 1, nvar)
          else if (mkind == 8) then
             inquire (iolength = nr) (1_I8_, ivar = 1, nvar)
          else
             inquire (iolength = nr) (1, ivar = 1, nvar)
          end if
       else
          inquire (iolength = nr) (1, ivar = 1, nvar)
       end if
    else                      ! character
       inquire (iolength = nr) ('1', ivar = 1, nvar)
    end if

  end function recordLen_I


  !+
  ! Read in 2-D array data in (x,y)=(col,row) format, from file
  !  Features automatic header skipping, first columns skipping, and automatic 
  !  detection of data sequence size 
  !-
  subroutine row_read_o2R(iu, nx, ny, dat, nskip, nyread)

    integer,  intent(in)  :: iu       ! file unit index
    integer,  intent(in)  :: nx       ! # of columns to be read
    integer,  intent(in)  :: ny       ! # of rows    to be read (0=unknown)
    real(R_), intent(out) :: dat(:,:) ! data read: Read data will be dat(nx, nyread)
    integer,  intent(in),  optional :: nskip  ! # of columns to skip (0 or larger)
    integer,  intent(out), optional :: nyread ! # of rows actualy read (if nyread < ny, error!)
    real(R_), allocatable :: wdum(:)
    integer :: ios, iy, nyqry, ny1, ns, ix, is

    ! Setup
    ns = 0
    if (present(nskip)) then
       ns = nskip
       allocate (wdum(ns))
    end if
    if (ny >= 1) then
       nyqry = ny      ! # of query data rows
    else
       nyqry = size(dat,2)
    end if
    ny1 = 0

    ! The first row
    do
       if (ns <= 0) then
          read (iu, *, iostat=ios) (dat(ix, 1), ix = 1, nx)
       else
          read (iu, *, iostat=ios) (wdum(is), is = 1, ns), (dat(ix, 1), ix = 1, nx)
       end if
       if (ios <= 0) exit
    end do
    if (ios < 0) then
       call err_read(ios, iu, 'row_read_o2R : No data found.')
       if (present(nyread)) nyread = ny1
       if (allocated(wdum)) deallocate(wdum)
       return
    end if
    ny1 = 1

    ! Read in subsequent rows
    if (ns <= 0) then
       do iy = 2, nyqry
          read (iu, *, iostat=ios) (dat(ix, iy), ix = 1, nx)
          if (ios /= 0) exit
          ny1 = ny1 + 1
       end do
    else
       do iy = 2, nyqry
          read (iu, *, iostat=ios) (wdum(is), is = 1, ns), (dat(ix, iy), ix = 1, nx)
          if (ios /= 0) exit
          ny1 = ny1 + 1
       end do
    end if
    if (ny1 < ny) call err_read(ios, iu, 'row_read_o2R : Too few data.')
    if (present(nyread)) nyread = ny1
    if (allocated(wdum)) deallocate(wdum)

  end subroutine row_read_o2R


  !+
  ! Read in 3-D array data in (x,y,z)=(col,row,lay) format, from file
  !  Features automatic header skipping, first columns skipping, and automatic 
  !  detection of data sequence size 
  !-
  subroutine row_read_o3R(iu, nx, ny, nz, dat, nskip, nzread)

    integer,  intent(in)  :: iu       ! file unit index
    integer,  intent(in)  :: nx       ! # of columns to be read
    integer,  intent(in)  :: ny       ! # of rows    to be read
    integer,  intent(in)  :: nz       ! # of layers  to be read (0=unknown)
    real(R_), intent(out) :: dat(:,:,:) ! data read: Read data will be dat(nx, ny, nzread)
    integer,  intent(in),  optional :: nskip  ! # of columns to skip (0 or larger)
    integer,  intent(out), optional :: nzread ! # of layers actualy read (if nzread < nz, error!)
    real(R_), allocatable :: wdum(:)
    integer :: ios, iy, iz, nzqry, nz1, iy_s, ns, ix, is

    ! Setup
    ns = 0
    if (present(nskip)) then
       ns = nskip
       allocate (wdum(ns))
    end if
    if (nz >= 1) then
       nzqry = nz      ! # of query data rows
    else
       nzqry = size(dat,2)
    end if
    nz1 = 0

    ! The first row
    do
       if (ns <= 0) then
          read (iu, *, iostat=ios) (dat(ix, 1, 1), ix = 1, nx)
       else
          read (iu, *, iostat=ios) (wdum(is), is = 1, ns), (dat(ix, 1, 1), ix = 1, nx)
       end if
       if (ios <= 0) exit
    end do
    if (ios < 0) then
       call err_read(ios, iu, 'row_read_o3R : No data found.')
       if (present(nzread)) nzread = nz1
       if (allocated(wdum)) deallocate(wdum)
       return
    end if

    ! Read in subsequent rows
    nz1 = 0
    iy_s = 2
    loop_z: do iz = 1, nzqry
       if (iz >= 2) iy_s = 1
       loop_y: do iy = iy_s, ny
          if (ns <= 0) then
             read (iu, *, iostat=ios) (dat(ix, iy, iz), ix = 1, nx)
          else
             read (iu, *, iostat=ios) (wdum(is), is = 1, ns), (dat(ix, iy, iz), ix = 1, nx)
          end if
          if (ios /= 0) exit loop_z
       end do loop_y
       nz1 = nz1 + 1
    end do loop_z
    if (nz1 < nz) call err_read(ios, iu, 'row_read_o3R : Too few data.')
    if (present(nzread)) nzread = nz1
    if (allocated(wdum)) deallocate(wdum)

  end subroutine row_read_o3R


  !+
  ! Write out 2-D array data in (x,y)=(col,row) format, to a file
  !  Features optional indexing and formatting
  !-
  subroutine row_write_i2R(iu, dat, nx, ny, midx, ncol, fmt)

    integer,  intent(in) :: iu         ! file unit index
    real(R_), intent(in) :: dat(:,:) ! data to be written: dat(nx, ny, nz)
    integer,  intent(in), optional :: nx    ! # of columns to be written
    integer,  intent(in), optional :: ny    ! # of rows    to be written
    integer,  intent(in), optional :: midx  ! method for indexing (0:no, 1:yes)
    integer,  intent(in), optional :: ncol  ! # of data columns in a single formatted line
    character(*), intent(in), optional :: fmt  ! format of writing a single datum value (e.g., 'es12.4')
    character(80) :: fmtstr
    integer :: nx1, ny1, midx1, ios, iy, ix

    ! Setup
    nx1 = size(dat, 1)
    ny1 = size(dat, 2)
    if (present(nx)) nx1 = nx
    if (present(ny)) ny1 = ny
    midx1 = 0
    if (present(midx)) midx1 = midx
    if (present(ncol) .and. present(fmt)) then
       if (midx1 <= 0) then
          fmtstr = trim(num2str_AN(ncol))//trim(fmt)
       else
          fmtstr = 'i6, '//trim(num2str_AN(nx1))//trim(fmt) ! ncol is neglected
       end if
       fmtstr = '('//trim(fmtstr)//')'
    else
       fmtstr = ' '
    end if

    ! Read in subsequent rows
    do iy = 1, ny1
       if (midx1 <= 0) then
          if (fmtstr == ' ') write (iu, *,      iostat=ios) (dat(ix, iy), ix = 1, nx1)
          if (fmtstr /= ' ') write (iu, fmtstr, iostat=ios) (dat(ix, iy), ix = 1, nx1)
       else
          if (fmtstr == ' ') write (iu, *,      iostat=ios) iy, (dat(ix, iy), ix = 1, nx1)
          if (fmtstr /= ' ') write (iu, fmtstr, iostat=ios) iy, (dat(ix, iy), ix = 1, nx1)
       end if
       if (ios /= 0) exit
    end do
    if (ios /= 0) call err_write(ios, iu, 'row_write_i2R :')

  end subroutine row_write_i2R


  !+
  ! Write out 3-D array data in (x,y,z)=(col,row,lay) format, to a file
  !  Features optional indexing and formatting
  !-
  subroutine row_write_i3R(iu, dat, nx, ny, nz, midx, ncol, fmt)

    integer,  intent(in) :: iu         ! file unit index
    real(R_), intent(in) :: dat(:,:,:) ! data to be written: dat(nx, ny, nz)
    integer,  intent(in), optional :: nx    ! # of columns to be written
    integer,  intent(in), optional :: ny    ! # of rows    to be written
    integer,  intent(in), optional :: nz    ! # of layers  to be written
    integer,  intent(in), optional :: midx  ! method for indexing (0:no, 1:yes)
    integer,  intent(in), optional :: ncol  ! # of data columns in a single formatted line
    character(*), intent(in), optional :: fmt  ! format of writing a single datum value (e.g., 'es12.4')
    character(80) :: fmtstr
    integer :: nx1, ny1, nz1, midx1, ios, iy, iz, ix

    ! Setup
    nx1 = size(dat, 1)
    ny1 = size(dat, 2)
    nz1 = size(dat, 3)
    if (present(nx)) nx1 = nx
    if (present(ny)) ny1 = ny
    if (present(nz)) nz1 = nz
    midx1 = 0
    if (present(midx)) midx1 = midx
    if (present(ncol) .and. present(fmt)) then
       if (midx1 <= 0) then
          fmtstr = trim(num2str_AN(ncol))//trim(fmt)
       else
          fmtstr = '2(1x, i6), '//trim(num2str_AN(nx1))//trim(fmt) ! ncol is neglected
       end if
       fmtstr = '('//trim(fmtstr)//')'
    else
       fmtstr = ' '
    end if

    ! Read in subsequent rows
    loop_z: do iz = 1, nz1
       do iy = 1, ny1
          if (midx1 <= 0) then
             if (fmtstr == ' ') write (iu, *,      iostat=ios) (dat(ix, iy, iz), ix = 1, nx1)
             if (fmtstr /= ' ') write (iu, fmtstr, iostat=ios) (dat(ix, iy, iz), ix = 1, nx1)
          else
             if (fmtstr == ' ') write (iu, *,      iostat=ios) iy, iz, (dat(ix, iy, iz), ix = 1, nx1)
             if (fmtstr /= ' ') write (iu, fmtstr, iostat=ios) iy, iz, (dat(ix, iy, iz), ix = 1, nx1)
          end if
          if (ios /= 0) exit loop_z
       end do
    end do loop_z
    if (ios /= 0) call err_write(ios, iu, 'row_write_i3R :')

  end subroutine row_write_i3R


  !+
  ! Skip lines in a sequential file by N lines
  !-
  subroutine skipLines(iu, nskip)

    integer, intent(in) :: iu    ! file unit index
    integer, intent(in) :: nskip ! # of lines to be skipped
    integer :: ios, i

    do i = 1, nskip
       read (iu, *, iostat=ios)
       call err_read(ios, iu)
    end do

  end subroutine skipLines


  !+
  ! Skip reading a text file until a user-specified sign (keyword) is found
  !-
  subroutine skipToSign(iu, str, ios)

    integer,      intent(in)  :: iu  ! file unit index
    character(*), intent(in)  :: str ! keywords
    integer,      intent(out) :: ios ! I/O status (0=normal)
    character(len_trim(str)) :: bufstr

    do 
       read (iu, '(a)', iostat=ios) bufstr
       if (ios /= 0) exit
       if (bufstr == trim(str)) exit
    end do

  end subroutine skipToSign


  !+
  ! Skip reading a text file until one of user-specified signs (keywords) is found
  !-
  subroutine skipToSigns(iu, strs, ios)

    integer,      intent(in)  :: iu      ! file unit index
    character(*), intent(in)  :: strs(:) ! keywords
    integer,      intent(out) :: ios     ! I/O status (0=normal)
    character(80) :: bufstr
    integer :: is, nc

    loop_skip: do 
       read (iu, '(a)', iostat=ios) bufstr
       if (ios /= 0) exit loop_skip
       do is = 1, size(strs)
          nc = len_trim(strs(is))
          if (bufstr(1:nc) == strs(is)(1:nc)) exit loop_skip
       end do
    end do loop_skip

  end subroutine skipToSigns

end module hparx_file
