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
! Library of utitilies that use nonstandard Fortran 90/95 procedures and functions
!
! -Note this module does not comform the Fortran 90/95 standard.
! -Fortran 2003 features are used in the current version.
!-
module hparx_nons

  use hparx_base, only : err_issue
  implicit none
  private

  public :: getCmdArgs ! Get command line arguments

contains

  !+
  ! Get command line arguments
  !-
  subroutine getCmdArgs(narg, argv, argmsg)

    integer, intent(inout) :: narg ! (in) least # of arguments, (out) actual # of arguments
    character(*), intent(out) :: argv(:) ! command arguments
    character(*), intent(in), optional :: argmsg ! error message to be printed
    integer   :: na, i

    ! Get # of arguments
    na = command_argument_count() ! Fortran 2003 and later

    ! Check # of arguments
    if (na < narg) then ! too few arguments
       if (present(argmsg)) then
          call err_issue(2, argmsg)
       else
          call err_issue(2, 'Too few command-line arguments.')
       end if
    endif
    narg = na

    ! Get arguments
    do i = 1, narg
       call get_command_argument(i, argv(i)) ! Fortran 2003 and later
    end do

  end subroutine getCmdArgs


  !+
  ! Get command line arguments
  !
  ! On the use of legacy getarg() and iargc()
  !  A procedure getarg and a function iargc are quite common but non-standard.
  !  Compilers may assume they are intrinsic or externally defined.
  !-
  !subroutine getCmdArgs(narg, argv, argmsg)
  !
  !  integer, intent(inout) :: narg ! (in) least # of arguments, (out) actual # of arguments
  !  character(*), intent(out) :: argv(:) ! command arguments
  !  character(*), intent(in), optional :: argmsg ! error message to be printed
  !  integer   :: na, i, id
  !
  !// For the PGI compiler, the following line should be activated:
  !  integer,  external :: iargc
  !// For other compilers (GNU, G95, Intel), the above should be deactivated.
  ! It may be useful to use sed to preprocess the code: (untested sed command)
  !   $ sed '/external.*iargc/d' < hparx_nons.F90 > hparx_nons.f90
  !
  ! Get # of arguments
  !  na = iargc() ! non-standard Fortran
  !  call getarg(0, argv(1)) ! non-standard Fortran
  !  if(argv(1) == "") then
  !     na = na - 1 ! for unusual compilers
  !     id = 1
  !  else
  !     id = 0 ! for usual compilers
  !  end if
  !
  ! Check # of arguments
  !  if (na < narg) then ! too few arguments
  !     if (present(argmsg)) then
  !        call err_issue(2, argmsg)
  !     else
  !        call err_issue(2, 'Too few command-line arguments.')
  !     end if
  !  endif
  !  narg = na
  !
  ! Get arguments
  !  do i = 1, narg
  !     call getarg(i + id, argv(i)) ! non-standard Fortran
  !  end do
  !
  !end subroutine getCmdArgs

end module hparx_nons
