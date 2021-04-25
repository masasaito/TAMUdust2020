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
! HPARX - a nice library of Fortran codes
!
! High Performance code library for the Atmospheric Sciences,
!  developed by the Radiation and climate physiX group
!  2020/04/29 MS this library includes only codes needed for
!  TAMUdust2020 database
!-
module hparx 

  use hparx_base
  use hparx_file
  use hparx_lina
  use hparx_math
  use hparx_vecmat
  use hparx_nons !* Includes procedures that do not conform Fortran standard *

  implicit none
  public  ! all in sub-modules are public
  !// This module is just a wrapper to submodules.

contains

  !+
  ! Print a description of this library
  !-
  subroutine hparx_help() 

    write (*,*) 'hparx, a nice library of Fortran codes'

  end subroutine hparx_help

end module hparx
