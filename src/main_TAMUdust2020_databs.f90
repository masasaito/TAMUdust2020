!+
! TAMUdust2020 
! 
!   Create a single-scattering property databse of user-defined aerosol particles 
!   based on the TAMUdust2020 data kernels
!
!   Update: 
!     Masa Saito  10/28/2019 Created an initial version
!     Masa Saito  07/31/2020 Released as v0.9
!     Masa Saito  09/27/2020 Released as v1.0 (Bugfix, implement the kernel method)
!     Masa Saito  04/21/2021 Multiple size definition incorporated.
!     Masa Saito  06/21/2021 Module structure incorporated.
!     Masa Saito  10/11/2021 Improved multi-dimensional linear interpolation.
! 
!-
program tamudust2020databs

  use tamudust2020, only : TAMUdust2020databs__init, TAMUdust2020__core, TAMUdust2020__final
  implicit none

  ! TAMUdust2020 Core Process
  call TAMUdust2020databs__init()
  call TAMUdust2020__core()
  call TAMUdust2020__final()
  stop

end program tamudust2020databs

