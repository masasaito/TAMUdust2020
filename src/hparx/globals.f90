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
! Global constants
!-
module globals 

  implicit none
  public            ! all in this module are public!

  ! Programmer-dependent kinds of real/complex variables
  !// Note: The R_ is default and RD_ should be used only when needed.
  !    One's choice of precision for R_ and RD_ may depend on required accuracy and used memory size.
  !integer,   parameter :: R_  = selected_real_kind(6)  ! default precision
  integer,   parameter :: R_  = selected_real_kind(13) ! default precision
  integer,   parameter :: RD_ = selected_real_kind(13) ! higher  precision

  ! System-specific kinds of variables (no freedom for programmers)
  integer,   parameter :: R4_ = selected_real_kind(6)  ! 4-byte real
  integer,   parameter :: R8_ = selected_real_kind(13) ! 8-byte real
  integer,   parameter :: I4_ = selected_int_kind(9)   ! 4-byte integer
  integer,   parameter :: I8_ = selected_int_kind(18)  ! 8-byte integer

  ! Numerical parameters for real variables
  real(R_),  parameter :: REPS_   = epsilon(1.0_R_)    ! 1.0e-7, typically
  real(R_),  parameter :: RTINY_  = tiny(1.0_R_)       ! 1.0e-37
  real(R_),  parameter :: RHUGE_  = huge(1.0_R_)       ! 1.0e+37
  real(R_),  parameter :: RSML_   = RTINY_  * 3.0_R_   ! 1.0e-37
  real(R_),  parameter :: RLRG_   = RHUGE_  / 3.0_R_   ! 1.0e+37
  real(RD_), parameter :: RDEPS_  = epsilon(1.0_RD_)   ! 1.0e-15, typically
  real(RD_), parameter :: RDTINY_ = tiny(1.0_RD_)      ! 1.0e-307
  real(RD_), parameter :: RDHUGE_ = huge(1.0_RD_)      ! 1.0e+307
  real(RD_), parameter :: RDSML_  = RDTINY_ * 3.0_RD_  ! 1.0e-307
  real(RD_), parameter :: RDLRG_  = RDHUGE_ / 3.0_RD_  ! 1.0e+307

  ! Numerical parameters for integer variables
  integer,   parameter :: IMAX_ =  2147483647          ! max integer of default kind
  integer,   parameter :: IMIN_ = -2147483647 - 1      ! min integer of default kind

  ! Mathematical constants
  real(R_),  parameter :: PI_     = 3.141592653589793238462643383279502884197_R_   ! pi
  real(R_),  parameter :: PIH_    = 1.570796326794896619231321691639751442098_R_   ! pi/2
  real(R_),  parameter :: PI2_    = 6.283185307179586476925286766559005768394_R_   ! pi*2
  real(RD_), parameter :: DPI_    = 3.141592653589793238462643383279502884197_RD_  ! pi
  real(RD_), parameter :: DPIH_   = 1.570796326794896619231321691639751442098_RD_  ! pi/2
  real(RD_), parameter :: DPI2_   = 6.283185307179586476925286766559005768394_RD_  ! pi*2
  real(R_),  parameter :: SQRT2_  = 1.414213562373095048801688724209698078569_R_   ! sqrt(2)
  real(R_),  parameter :: SQRT3_  = 1.732050807568877293527446341505872366942_R_   ! sqrt(3)
  real(R_),  parameter :: SQRT2A_ = 0.7071067811865475244008443621048490392848_R_  ! 1/sqrt(2)
  real(R_),  parameter :: SQRT3A_ = 0.5773502691896257645091487805019574556476_R_  ! 1/sqrt(3)
  real(R_),  parameter :: DTOR_   = PI_ / 180.0_R_     ! degree-to-radian factor
  real(R_),  parameter :: FRAC13_ = 1.0_R_ / 3.0_R_    ! 1/3
  real(R_),  parameter :: FRAC23_ = 2.0_R_ / 3.0_R_    ! 2/3
  real(R_),  parameter :: FRAC43_ = 4.0_R_ / 3.0_R_    ! 4/3
  real(R_),  parameter :: FRAC53_ = 5.0_R_ / 3.0_R_    ! 5/3
  real(R_),  parameter :: FRAC16_ = 1.0_R_ / 6.0_R_    ! 1/6

  ! Physical constants (SI unit)
  real(R_),  parameter :: CSPEED_ = 2.99792458e+8_R_     ! speed of light (m/s)
  real(R_),  parameter :: PLANCK_ = 6.6260689633e-34_R_  ! Planck's constant (J*s)
  real(R_),  parameter :: BOLTZ_  = 1.380650424e-23_R_   ! Boltzmann's constant (J/K)
  real(R_),  parameter :: STEBOL_ = 5.67040040e-8_R_     ! Stefan-Boltzmann's constant (W/m^2/K^4)
  real(R_),  parameter :: AVOGAD_ = 6.0221417930e+23_R_  ! Avogadro number (/mol)
  real(R_),  parameter :: GASCON_ = 8.31447215_R_        ! universal gas constant (J/mol/K)
  real(R_),  parameter :: VOLMOL_ = 2.241399639e-2_R_    ! molar volume (m^3/mol) of ideal gas (0 deg.C, 1 atm)
  real(R_),  parameter :: EGRAVE_ = 9.7803267715_R_      ! earth equatorial gravity (m/s^2)
  real(R_),  parameter :: EGRAV0_ = 9.80665_R_           ! earth standard   gravity (m/s^2)
  real(R_),  parameter :: WGTVAP_ = 18.016e-3_R_         ! molar mass of water vapor    (H2O) (kg/mol)
  real(R_),  parameter :: WGTCO2_ = 44.01e-3_R_          ! molar mass of carbon dioxide (CO2) (kg/mol)
  real(R_),  parameter :: WGTOZN_ = 48.00e-3_R_          ! molar mass of ozone          (O3)  (kg/mol)
  real(R_),  parameter :: WGTDRY_ = 28.96454e-3_R_       ! molar mass of dry air (kg/mol)
  real(R_),  parameter :: RHOWAT0_ = 0.99984e+3_R_       ! density (kg/m^3) of liquid water (0 deg.C)
  real(R_),  parameter :: RHOICE0_ = 0.917e+3_R_         ! density (kg/m^3) of water ice (0 deg.C)

contains

  !+
  ! Print values of some global constants
  !-
  subroutine globals__print() 

    write (*,*) R_, RD_, R4_, R8_
    write (*,*) I4_, I8_
    write (*,*) REPS_, RTINY_, RHUGE_
    write (*,*) IMIN_, IMAX_

  end subroutine globals__print

end module globals
