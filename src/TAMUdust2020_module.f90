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
module tamudust2020

  use globals
  use hparx, only : open_dir, gridIdx_loc, close_all, open_seq, getCmdArgs, err_read, &
                    intpTab_linear_1R
  implicit none
  private
  public :: TAMUdust2020databs__init
  public :: TAMUdust2020single__init
  public :: TAMUdust2020__core
  public :: TAMUdust2020__final


  ! TAMUdust2020 binary kernel table
  integer, save      :: Tab_nsiz   = 169      ! # of size bin
  integer, parameter :: Tab_nang   = 498      ! # of scattering angle bin
  integer, parameter :: Tab_npol   = 6        ! # of polarization elements
  integer, parameter :: Tab_nmet   = 9        ! # of meta information (spara, Aproj, Vol, Cext, Csca, g*Csca, RefRe, RefImi, Sphr) 
  integer, save      :: Tab_nsph   = 6        ! # of sphericity bin
  integer, save      :: Tab_nReIdx = 8        ! # of real part of dust refractive index bin
  integer, save      :: Tab_nImIdx = 31       ! # of imaginary part of dust reflactive index bin
  integer, save      :: Tab_nref              ! # of refractive index bin
  integer, save      :: Tab_ReIdxSoft_ID = 1  ! ID for optically soft refractive index (=1.1)
  integer, parameter :: Tab_nfile  = 2        ! binary pmat and isca files
  real(R4_), allocatable :: TAMUdust2020STab(:,:,:,:) ! isca.dat kernel table
  real(R4_), allocatable :: TAMUdust2020PTab(:,:,:,:) ! Phase matrix dat kernel table 
  real(R4_), allocatable :: Tab_siz(:), Tab_ReIdx(:), Tab_ImIdx(:), Tab_sph(:) 

  ! TAMUdust2020 ASCII user-defined aerosol scattering property database
  character(len=256)  :: Dust_outpath       ! output dust database path
  character(len=10)   :: Dust_siztyp        ! size type definition
  character(len=9)    :: Dust_dattyp        ! data type (DataTable, SingleCal, ErrorEval)
  integer             :: Dust_nfile         ! # of files
  integer, parameter  :: Dust_nmet  = 7     ! (wav, Dmax, Vol, Aproj, Qext, SSA, g)
  integer, parameter  :: Dust_npol  = 6     ! (P11, P12, P22, P33, P43, P44)
  integer, save       :: Dust_nsiz          ! # of size bin 
  integer, save       :: Dust_nang          ! # of scattering angle
  integer, save       :: Dust_nwav          ! # of wavelength
  real(R_), allocatable :: Dust_siz(:)      ! (nsiz) maximum dimension
  real(R_), allocatable :: Dust_refidx(:,:) ! (nref,2) refractive index (Re,Im)
  real(R_), allocatable :: Dust_ang(:)      ! (nang) scattering angle
  real(R_), allocatable :: Dust_wav(:)      ! (nwav) Wavelength
  real(R_), save        :: Dust_sph         ! Sphericity

  ! Numerical parameter
  real(R4_), parameter :: PI4_ = 3.14159265358_R4_

contains


   !+
   ! TAMUdust2020: Initialization for database creation
   !-
   subroutine TAMUdust2020databs__init()

      integer, parameter :: IUI = 111, IUO = 112
      integer            :: isiz, iwav, iang, iref, ire, iim, isph
      integer            :: irec, irecls, ireclp, ios, ioss, iosp, ivar, narg
      character(len=256) :: argv(10) ! argument to be red
      character(len=256) :: argmsg   ! argment message
      character(len=256), save :: sizfile           ! (nsiz) maximum dimension file
      character(len=256), save :: wavfile           ! (nwav) wavelength file
      character(len=256), save :: refidxfile        ! (nwav,2) refractive index file
      character(len=256), save :: angfile           ! (nang) scattering angle
      character(len=256), save :: outpath           ! output path of horizontally oriented phase matrix
      character(len=10),  save :: siztyp = 'MaximumDim' ! size type definition
      character(len=256)       :: TAMUdust2020path  ! TAMUdust2020 binary table path
      integer,  save           :: nsiz, nwav        ! # of size/wavelength bins
      integer, save            :: nsiztab, nrefretab, nrefimtab, nsphtab
      real(R_), save           :: spher
      namelist /TAMUdust2020/TAMUdust2020path, outpath, sizfile, wavfile, refidxfile, angfile, &
                             siztyp, nsiz, nwav, spher, nsiztab, nrefretab, nrefimtab, nsphtab            

      argmsg = ' Usage: tamudust2020databs nmlst > output '  
      narg = 1
      call getCmdArgs(narg, argv, argmsg)
      open(IUI, file = trim(adjustl(argv(1))), status = 'old')
      read(IUI, nml=TAMUdust2020, iostat=ios)
      call err_read(ios, IUI, 'Main: Error in reading the namelist data.')
      close(IUI)


      ! Pre-setting
      Dust_dattyp = 'DataTable'
      Dust_nfile  = 8


      ! Initialization
      Dust_outpath = outpath
      Dust_siztyp  = siztyp
      Dust_nsiz    = nsiz
      Dust_nwav    = nwav
      Dust_nang    = Tab_nang  ! This will be flexible 
      Dust_sph     = spher
      Tab_nsiz     = nsiztab
      Tab_nReIdx   = nrefretab
      Tab_nImIdx   = nrefimtab
      Tab_nsph     = nsphtab
      Tab_nref     = Tab_nReIdx * Tab_nImIdx
      allocate( Dust_siz(Dust_nsiz), Dust_refidx(Dust_nwav,2), Dust_ang(Dust_nang), Dust_wav(Dust_nwav))
      allocate( TAMUdust2020STab(Tab_nmet,         Tab_nsph,Tab_nsiz,Tab_nref) ) ! isca.dat
      allocate( TAMUdust2020PTab(Tab_nang*Tab_npol,Tab_nsph,Tab_nsiz,Tab_nref) ) ! Phase matrix dat 
      allocate( Tab_siz(Tab_nsiz), Tab_ReIdx(Tab_nReIdx), Tab_ImIdx(Tab_nImIdx), Tab_sph(Tab_nsph) ) 


      ! Read files
      open(IUI, file = trim(sizfile), status = 'old') ! size file
      do isiz = 1, Dust_nsiz
         read(IUI,*) Dust_siz(isiz)
      end do
      close(IUI)
      open(IUI, file = trim(wavfile), status = 'old') ! wavelength file
      do iwav = 1, Dust_nwav
         read(IUI,*) Dust_wav(iwav)
      end do
      close(IUI)
      open(IUI, file = trim(refidxfile), status = 'old') ! refidx file
      do iwav = 1, Dust_nwav
         read(IUI,*) Dust_refidx(iwav,1:2)
      end do
      close(IUI)
      open(IUI, file = trim(angfile), status = 'old') ! refidx file
      do iang = 1, Dust_nang
         read(IUI,*) Dust_ang(iang)
      end do
      close(IUI)
 

      ! Inquire and load database
      inquire (iolength=irecls) (1.0_R4_, ivar=1, Tab_nmet)
      inquire (iolength=ireclp) (1.0_R4_, ivar=1, Tab_nang*Tab_npol)
      call open_dir(IUI, trim(TAMUdust2020path)//'/TAMUdust2020isca.bin', irecls, 'old')
      call open_dir(IUO, trim(TAMUdust2020path)//'/TAMUdust2020PMat.bin', ireclp, 'old')
      irec = 1
      do isiz = 1, Tab_nsiz
         do iref = 1, Tab_nref
            do isph = 1, Tab_nsph
               read(IUI, rec=irec, iostat=ioss) TAMUdust2020STab(:,isph,isiz,iref) 
               read(IUO, rec=irec, iostat=iosp) TAMUdust2020PTab(:,isph,isiz,iref) 
               if (ioss /= 0 .or. iosp /= 0) stop 'ERROR: database is something wrong'
               irec = irec + 1
            end do
         end do
      end do
      close(IUI)
      close(IUO)
    
      ! Size type definition
      if (Dust_siztyp /= 'MaximumDim') then ! not Dmax
         if      (Dust_siztyp == 'ProjecEqvS') then ! Dare = (4*A/pi)^1/2
            TAMUdust2020STab(1,:,:,:) = (TAMUdust2020STab(2,:,:,:)*4.0_R4_/PI4_)**0.5_R4_
         else if (Dust_siztyp == 'VolumeEqvS') then ! Dvol = (6*V/pi)^1/3
            TAMUdust2020STab(1,:,:,:) = (TAMUdust2020STab(3,:,:,:)*6.0_R4_/PI4_)**0.3333333333_R4_
         else
            stop 'siztyp should be either ("MaximumDim", "ProjecEqvS", "VolumeEqvS")'
         end if
      endif


      ! Get meta information
      Tab_ReIdxSoft_ID = 1
      Tab_siz(1:Tab_nsiz) =  TAMUdust2020STab(1,1,1:Tab_nsiz,1)
      do ire = 1, Tab_nReIdx
         Tab_ReIdx(ire) =  TAMUdust2020STab(7,1,1,1+(ire-1)*Tab_nImIdx)
         if (Tab_ReIdx(ire) < 1.1005_R4_) Tab_ReIdxSoft_ID = ire
         do iim = 1, Tab_nImIdx
            Tab_ImIdx(iim) =  TAMUdust2020STab(8,1,1,iim)
         end do
      end do
      do isph = 1, Tab_nsph
         Tab_sph(isph) = TAMUdust2020STab(9,isph,1,1)
      end do 


      ! Warnings
      call TAMUdust2020__warning()


   end subroutine TAMUdust2020databs__init



   !+
   ! TAMUdust2020: Read namelist
   !-
   subroutine TAMUdust2020single__init()
     
      integer, parameter :: IUI = 111, IUO = 112
      integer            :: isiz, iang, iref, isph, ire, iim
      integer            :: irec, irecls, ireclp, ios, ioss, iosp, ivar, narg
      integer, save      :: nsiztab, nrefretab, nrefimtab, nsphtab
      character(len=256) :: argv(10)   ! argument to be red
      character(len=290) :: argmsg     ! argment message
      character(len=256), save :: angfile           ! (nang) scattering angle
      character(len=10),  save :: siztyp = 'MaximumDim' ! size type definition
      character(len=256)       :: TAMUdust2020path  ! TAMUdust2020 binary table path
      namelist /TAMUdust2020/TAMUdust2020path, angfile, siztyp, nsiztab, nrefretab, nrefimtab, nsphtab            

      ! Arguments
      narg = 7
      argmsg = 'Usage: tamudust2020single nmlst siz wav refRe refIm outpath \n\n' &
       &     //'  nmlst   : namelist file\n' &
       &     //'  siz     : particle size (um)\n' &
       &     //'  wav     : wavelength (um)\n' &
       &     //'  refRe   : real part of refractive index\n' &
       &     //'  refIm   : imaginary part of refractive index\n' &
       &     //'  sph     : sphericity\n' &
       &     //'  outpath : output path' 
      call getCmdArgs(narg, argv, argmsg)
      open(IUI, file = trim(adjustl(argv(1))), status = 'old')
      read(IUI, nml=TAMUdust2020, iostat=ios)
      call err_read(ios, IUI, 'Main: Error in reading the namelist data.')
      close(IUI)

  
      ! Pre-setting
      Dust_dattyp = 'SingleCal'
      Dust_siztyp  = siztyp
      Dust_nfile  = 2
      Dust_nsiz   = 1
      Dust_nwav   = 1
      Dust_nang    = Tab_nang  ! This will be flexible 
      allocate( Dust_siz(Dust_nsiz), Dust_refidx(Dust_nwav,2), Dust_ang(Dust_nang), Dust_wav(Dust_nwav) )
 

      ! Initialization
      read(argv(2), *)     Dust_siz(1)
      read(argv(3), *)     Dust_wav(1)
      read(argv(4), *)     Dust_refidx(1,1)
      read(argv(5), *)     Dust_refidx(1,2)
      read(argv(6), *)     Dust_sph
      read(argv(7), '(a)') Dust_outpath
      Tab_nsiz     = nsiztab
      Tab_nReIdx   = nrefretab
      Tab_nImIdx   = nrefimtab
      Tab_nsph     = nsphtab
      Tab_nref     = Tab_nReIdx * Tab_nImIdx
      allocate( TAMUdust2020STab(Tab_nmet,         Tab_nsph,Tab_nsiz,Tab_nref) ) ! isca.dat
      allocate( TAMUdust2020PTab(Tab_nang*Tab_npol,Tab_nsph,Tab_nsiz,Tab_nref) ) ! Phase matrix dat 
      allocate( Tab_siz(Tab_nsiz), Tab_ReIdx(Tab_nReIdx), Tab_ImIdx(Tab_nImIdx), Tab_sph(Tab_nsph) ) 


      ! Read files
      open(IUI, file = trim(angfile), status = 'old') ! refidx file
      do iang = 1, Dust_nang
         read(IUI,*) Dust_ang(iang)
      end do
      close(IUI)
 

      ! Inquire and load database
      inquire (iolength=irecls) (1.0_R4_, ivar=1, Tab_nmet)
      inquire (iolength=ireclp) (1.0_R4_, ivar=1, Tab_nang*Tab_npol)
      call open_dir(IUI, trim(TAMUdust2020path)//'/TAMUdust2020isca.bin', irecls, 'old')
      call open_dir(IUO, trim(TAMUdust2020path)//'/TAMUdust2020PMat.bin', ireclp, 'old')
      irec = 1
      do isiz = 1, Tab_nsiz
         do iref = 1, Tab_nref
            do isph = 1, Tab_nsph
               read(IUI, rec=irec, iostat=ioss) TAMUdust2020STab(:,isph,isiz,iref) 
               read(IUO, rec=irec, iostat=iosp) TAMUdust2020PTab(:,isph,isiz,iref) 
               if (ioss /= 0 .or. iosp /= 0) stop 'ERROR: database is something wrong'
               irec = irec + 1
            end do
         end do
      end do
      close(IUI)
      close(IUO)

      ! Size type definition
      if (Dust_siztyp /= 'MaximumDim') then ! not Dmax
         if      (Dust_siztyp == 'ProjecEqvS') then ! Dare = (4*A/pi)^1/2
            TAMUdust2020STab(1,:,:,:) = (TAMUdust2020STab(2,:,:,:)*4.0_R4_/PI4_)**0.5_R4_
         else if (Dust_siztyp == 'VolumeEqvS') then ! Dvol = (6*V/pi)^1/3
            TAMUdust2020STab(1,:,:,:) = (TAMUdust2020STab(3,:,:,:)*6.0_R4_/PI4_)**0.3333333333_R4_
         else
            stop 'siztyp should be either ("MaximumDim", "ProjecEqvS", "VolumeEqvS")'
         end if
      endif


      ! Get meta information
      Tab_ReIdxSoft_ID = 1
      Tab_siz(1:Tab_nsiz) = TAMUdust2020STab(1,1,1:Tab_nsiz,1)
      do ire = 1, Tab_nReIdx
         Tab_ReIdx(ire) =  TAMUdust2020STab(7,1,1,1+(ire-1)*Tab_nImIdx)
         if (Tab_ReIdx(ire) < 1.1005_R4_) Tab_ReIdxSoft_ID = ire
         do iim = 1, Tab_nImIdx
            Tab_ImIdx(iim) =  TAMUdust2020STab(8,1,1,iim)
         end do
      end do
      do isph = 1, Tab_nsph
         Tab_sph(isph) = TAMUdust2020STab(9,isph,1,1)
      end do 


      !Warnings
      call TAMUdust2020__warning()


   end subroutine TAMUdust2020single__init

 
   !+
   ! TAMUdust2020: Warnings
   ! ISSUE: Due to single precision, the mimimum refractive index is a bit
   ! larger than the minimum, providing warning when specify the minimum
   ! refractive index. (06/16/2021)
   !-
   subroutine TAMUdust2020__warning()
 
      ! Warnings
      if (maxval(real(Dust_siz))*2.0_R4_*PI4_/minval(real(Dust_wav)) > maxval(Tab_siz)) &
      print*, 'Warning: large sizes in the database may be based on extrapolation. Size parameter ', &
      maxval(real(Dust_siz))*2.0_R4_*PI4_/minval(real(Dust_wav)), ' is beyond the prescribed range.' 
      if (minval(Dust_refidx(:,1)) < minval(Tab_ReIdx(:))) &
      print*, 'Warning: Specified real part of refractive index', minval(Dust_refidx(:,1)),' is beyond the prescribed range. '
      if (maxval(Dust_refidx(:,1)) > maxval(Tab_ReIdx(:))) &
      print*, 'Warning: Specified real part of refractive index', maxval(Dust_refidx(:,1)),' is beyond the prescribed range. '
      if (minval(Dust_refidx(:,2)) < minval(Tab_ImIdx(:))) &
      print*, 'Warning: Specified imaginary part of refractive index', minval(Dust_refidx(:,2)),' is beyond the prescribed range. '
      if (maxval(Dust_refidx(:,2)) > maxval(Tab_ImIdx(:))) &
      print*, 'Warning: Specified imaginary part of refractive index', maxval(Dust_refidx(:,2)),' is beyond the prescribed range. '
      if (Dust_sph < minval(Tab_sph(:))  .or. Dust_sph > maxval(Tab_sph(:))) &
      print*, 'Warning: Specified sphericity ', Dust_sph, ' is beyond the prescribed range. '

   end subroutine TAMUdust2020__warning
  


   !+
   ! TAMUdust2020: Finalize the database
   !-
   subroutine TAMUdust2020__final()

      deallocate( Dust_siz, Dust_refidx, Dust_ang, Dust_wav )
      deallocate( TAMUdust2020STab, TAMUdust2020PTab, Tab_siz, Tab_ReIdx, Tab_ImIdx, Tab_sph ) 

   end subroutine TAMUdust2020__final


   !+
   ! TAMUdust2020: Core process based on user-defined configuration 
   !-
   subroutine TAMUdust2020__core()
     
      real(R_), parameter :: efac(Tab_nmet) = (/ 1.0_R_, 2.0_R_, 3.0_R_, 2.0_R_, & ! factor for linear interpol
                                         2.0_R_, 2.0_R_, 1.0_R_, 1.0_R_, 1.0_R_ /)
      integer  :: iwav, isiz, imet, ire, iim, isp
      real(R_) :: sre, sim, ssp 
      integer,  allocatable :: iuo(:)        ! I/O unit
      real(R_), allocatable :: wrk1s(:,:), wrk1p(:,:)      ! wrk for isca/p??.dat.    
      real(R_), allocatable :: wrk2s(:,:), wrk2p(:,:)      ! wrk for isca/p??.dat.    
      real(R_), allocatable :: isca(:,:),  pmat(:,:,:)     ! isca/p??.dat.    
      real(R_), allocatable :: spara(:)

      ! Set properties
      allocate( iuo(Dust_nfile), spara(Dust_nsiz) )
      allocate( wrk1s(Tab_nmet,Tab_nsiz),  wrk1p(Tab_nang*Tab_npol,Tab_nsiz)     )
      allocate( wrk2s(Tab_nmet,Dust_nsiz), wrk2p(Tab_nang*Tab_npol,Dust_nsiz)    )
      allocate( isca(Dust_nmet,Dust_nsiz), pmat(Dust_nang, Dust_npol, Dust_nsiz) )
      call TAMUdust2020_database_open(iuo)

      ! Interpolation
      do iwav = 1, Dust_nwav

         ! Pickup data based on a specified refractive index
         if      (Dust_refidx(iwav,1) > 1.1005_R_) then
            call gridIdx_loc(log(dble(Tab_ReIdx(Tab_ReIdxSoft_ID:))-1.0_R_), log(Dust_refidx(iwav,1)-1.0_R_), ire, sre, 1)
            ire = ire + Tab_ReIdxSoft_ID - 1
         else if (Dust_refidx(iwav,1) < 0.8995_R_) then
            call gridIdx_loc(log(1.0_R_-dble(Tab_ReIdx(:))), log(1.0_R_-Dust_refidx(iwav,1)), ire, sre, 1)
            sre = 1.0_R_ - sre ! 
         else
            call gridIdx_loc(dble(Tab_ReIdx(:)), dble(Dust_refidx(iwav,1)), ire, sre, 1)
         end if
         call gridIdx_loc(log(dble(Tab_ImIdx(:))), log(Dust_refidx(iwav,2)), iim, sim, 1)
         call gridIdx_loc(log(dble(Tab_sph(:))),   log(Dust_sph),            isp, ssp, 1)
         sre = max(0.0_R_, min(1.0_R_, sre))
         sim = max(0.0_R_, min(1.0_R_, sim))
         ssp = max(0.0_R_, min(1.0_R_, ssp))
         wrk1s(:,:) = (1.0_R_-sre)*(1.0_R_-sim)*(1.0_R_-ssp) * TAMUdust2020STab(:,isp  ,:,iim  +(ire-1)*Tab_nImIdx) &  
                    +         sre *(1.0_R_-sim)*(1.0_R_-ssp) * TAMUdust2020STab(:,isp  ,:,iim  + ire   *Tab_nImIdx) &
                    + (1.0_R_-sre)*        sim *(1.0_R_-ssp) * TAMUdust2020STab(:,isp  ,:,iim+1+(ire-1)*Tab_nImIdx) &
                    +         sre *        sim *(1.0_R_-ssp) * TAMUdust2020STab(:,isp  ,:,iim+1+ ire   *Tab_nImIdx) &
                    + (1.0_R_-sre)*(1.0_R_-sim)*        ssp  * TAMUdust2020STab(:,isp+1,:,iim  +(ire-1)*Tab_nImIdx) &  
                    +         sre *(1.0_R_-sim)*        ssp  * TAMUdust2020STab(:,isp+1,:,iim  + ire   *Tab_nImIdx) &
                    + (1.0_R_-sre)*        sim *        ssp  * TAMUdust2020STab(:,isp+1,:,iim+1+(ire-1)*Tab_nImIdx) &
                    +         sre *        sim *        ssp  * TAMUdust2020STab(:,isp+1,:,iim+1+ ire   *Tab_nImIdx) 
         wrk1p(:,:) = (1.0_R_-sre)*(1.0_R_-sim)*(1.0_R_-ssp) * TAMUdust2020PTab(:,isp  ,:,iim  +(ire-1)*Tab_nImIdx) &  
                    +         sre *(1.0_R_-sim)*(1.0_R_-ssp) * TAMUdust2020PTab(:,isp  ,:,iim  + ire   *Tab_nImIdx) &
                    + (1.0_R_-sre)*        sim *(1.0_R_-ssp) * TAMUdust2020PTab(:,isp  ,:,iim+1+(ire-1)*Tab_nImIdx) &
                    +         sre *        sim *(1.0_R_-ssp) * TAMUdust2020PTab(:,isp  ,:,iim+1+ ire   *Tab_nImIdx) &
                    + (1.0_R_-sre)*(1.0_R_-sim)*        ssp  * TAMUdust2020PTab(:,isp+1,:,iim  +(ire-1)*Tab_nImIdx) &  
                    +         sre *(1.0_R_-sim)*        ssp  * TAMUdust2020PTab(:,isp+1,:,iim  + ire   *Tab_nImIdx) &
                    + (1.0_R_-sre)*        sim *        ssp  * TAMUdust2020PTab(:,isp+1,:,iim+1+(ire-1)*Tab_nImIdx) &
                    +         sre *        sim *        ssp  * TAMUdust2020PTab(:,isp+1,:,iim+1+ ire   *Tab_nImIdx) 

         ! Particle size parameter determination
         spara(:) = 2.0_R_ * PI_ * Dust_siz(:) / Dust_wav(iwav)
         !spara(:) = PI_ * Dust_siz(:) / Dust_wav(iwav) ! 09/14/2020 BUGFIXED
         !if (minval(spara(:)) >= minval(Tab_siz(:))) then ! 04/21/2021 multiple size def. incorporated.
         if (minval(spara(:)) >= minval(wrk1s(1,:))) then
            isiz = 1
         else

            ! Single-scattering properties: Rayleigh scattering approximation
            loop_minsiz: do isiz = 1, Tab_nsiz  
               !if (spara(isiz) >= minval(Tab_siz(:))) exit loop_minsiz 
               ! 04/21/2021 multiple size def. incorporated.
               if (spara(isiz) >= minval(wrk1s(1,:))) exit loop_minsiz
            end do loop_minsiz          
            call TAMUdust2020_RayleighScat(spara(1:isiz-1), wrk1s(:,1), wrk2s(:,1:isiz-1), wrk2p(:,1:isiz-1))
         end if

         ! Single-scattering properties: multi-dimensional interpolation from the table
         do imet = 1, Tab_nmet
            !wrk2s(imet,isiz:) = intpTab_Akima_1R(log(Tab_siz(:)), wrk1s(imet,:), log(spara(isiz:)), 1)  
            ! 04/21/2021 multiple size def. incorporated.
            !wrk2s(imet,isiz:) = intpTab_Akima_1R(log(wrk1s(1,:)), wrk1s(imet,:), log(spara(isiz:)), 1)  
            ! 10/11/2021 more reliable interpolation
            wrk2s(imet,isiz:) = intpTab_linear_1R((wrk1s(1,:))**efac(imet), wrk1s(imet,:), (spara(isiz:))**efac(imet), 1)  
         end do   
         do imet = 1, Tab_nang*Tab_npol
            !wrk2p(imet,isiz:) = intpTab_Akima_1R(log(Tab_siz(:)), wrk1p(imet,:), log(spara(isiz:)), 1)  
            ! 04/21/2021 multiple size def. incorporated.
            !wrk2p(imet,isiz:) = intpTab_Akima_1R(log(wrk1s(1,:)), wrk1p(imet,:), log(spara(isiz:)), 1)  
            ! 10/11/2021 more reliable interpolation.
            wrk2p(imet,isiz:) = intpTab_linear_1R((wrk1s(1,:))**2.0_R_, wrk1p(imet,:), (spara(isiz:))**2.0_R_, 1)  
         end do    

         ! write 
         call TAMUdust2020_convert(wrk2s, wrk2p, iwav, isca, pmat)
         if (Dust_dattyp == 'DataTable') call TAMUdust2020_write_database(iuo, isca, pmat)
      end do
     
      ! Finalize the database
      if (Dust_dattyp == 'DataTable') then
         call TAMUdust2020_write_ReadMe(iuo)
      else if (Dust_dattyp == 'SingleCal') then
         call TAMUdust2020_write_single(iuo, isca, pmat)
      end if
      call close_all(iuo)
      deallocate( iuo, spara, wrk1s, wrk1p, wrk2s, wrk2p, isca, pmat)

   end subroutine TAMUdust2020__core


   !+
   ! TAMUdust2020: Rayleifgh scattering calculation based on V/A-equivalent sphere
   !-
   subroutine TAMUdust2020_RayleighScat(spara, wrk1s, wrk2s, wrk2p)
  
     real(R_), intent(in)  :: spara(:)   ! (nrsiz) Size parameter  
     real(R_), intent(in)  :: wrk1s(:)   ! (nmet) isca.dat from table  
                                         ! 1st column: Dmax, Aproj, Vol, Qext, SSA, g, RefRe, RefIm, Sphr    
     real(R_), intent(out) :: wrk2s(:,:) ! (nmet,nrsiz) isca.dat for database  
     real(R_), intent(out) :: wrk2p(:,:) ! (npol*nang,nrsiz) P??.dat for database  
     complex(R_) :: mref
     integer     :: isiz
     real(R_)    :: sfac, rva, mre, mim, qabs, qsca
 

     ! Size loop     
     do isiz = 1, size(spara)

        ! Parameters
        !sfac = spara(isiz)/Tab_siz(1)
        sfac = spara(isiz)/wrk1s(1) ! 04/21/2021 multiple size def. incorporated.
        mref = dcmplx(wrk1s(7), -wrk1s(8))
        mre  = real(((mref**2.0_R_-1.0_R_)/(mref**2.0_R_+2.0_R_))**2.0_R_)
        mim  = -aimag((mref**2.0_R_-1.0_R_)/(mref**2.0_R_+2.0_R_))
        rva  = 0.75_R_ * sfac * wrk1s(3) / wrk1s(2)              ! V/A-equivalent sphere radius (micron)
        qsca = 2.66666666_R_ * rva**4.0_R_ * mre
        qabs = 4.0_R_ * rva * mim

        ! isca.dat
        wrk2s(1,isiz)   = wrk1s(1)*sfac                       ! Dmax   (Dimensionless)
        wrk2s(2,isiz)   = wrk1s(2)*sfac**2.0_R_               ! Aproj  (Dimensionless)
        wrk2s(3,isiz)   = wrk1s(3)*sfac**3.0_R_               ! Vol    (Dimensionless)
        wrk2s(4,isiz)   = (qsca + qabs)*wrk1s(2)*sfac**2.0_R_ ! Cext   (Dimensionless)
        wrk2s(5,isiz)   = qsca*wrk1s(2)*sfac**2.0_R_          ! Csca   (Dimensionless)
        wrk2s(6,isiz)   = 0.0_R_                              ! Csca*g (Dimensionless)
        wrk2s(7:9,isiz) = wrk1s(7:9)                          ! RefRe, RefIm, Sphericity

        ! P??.dat
        wrk2p(           1:  Tab_nang,isiz) =  0.75_R_*(1.0_R_+dcos(Dust_ang(1:Tab_nang)*DTOR_)**2.0_R_) ! P11
        wrk2p(  Tab_nang+1:2*Tab_nang,isiz) = -0.75_R_*(1.0_R_-dcos(Dust_ang(1:Tab_nang)*DTOR_)**2.0_R_) ! P12
        wrk2p(2*Tab_nang+1:3*Tab_nang,isiz) =  wrk2p(           1:  Tab_nang,isiz)                       ! P22
        wrk2p(3*Tab_nang+1:4*Tab_nang,isiz) =  1.5_R_ * dcos(Dust_ang(1:Tab_nang)*DTOR_)                 ! P33
        wrk2p(4*Tab_nang+1:5*Tab_nang,isiz) =  0.0_R_                                                    ! P43
        wrk2p(5*Tab_nang+1:6*Tab_nang,isiz) =  wrk2p(3*Tab_nang+1:4*Tab_nang,isiz)                       ! P44
        wrk2p(           1:6*Tab_nang,isiz) =  wrk2p(           1:6*Tab_nang,isiz) * qsca*wrk1s(2)*sfac**2.0_R_ ! P?? --> Csca*P??
     end do
   
   end subroutine TAMUdust2020_RayleighScat


   !+ 
   ! TAMUdust2020: Format the database
   !-
   subroutine TAMUdust2020_convert(wrks, wrkp, iwav, isca, pmat)
 
      real(R_), intent(in)  :: wrks(:,:)     ! (nvar=9,nsiz)    isca.dat 
      real(R_), intent(in)  :: wrkp(:,:)     ! (nang*npol,nsiz) pmat.dat
      integer,  intent(in)  :: iwav
      real(R_), intent(out) :: isca(:,:)     ! (nvar=7,nsiz)    isca.dat
      real(R_), intent(out) :: pmat(:,:,:)   ! (nang,npol,nsiz) pmat.dat
      real(R_) :: sfac
      integer  :: ipol, iang
      !//NOTE wrks: first column (spara, Aproj, Vol, Cext, Csca, g*Csca, RefRe, RefIm, Sphr)
      ! for Aproj and Vol, scaling from size parameter to actual size is needed.

      ! convert optical properties
      !sfac      = Dust_wav(iwav) / PI_
      sfac      = Dust_wav(iwav) / (2.0_R_ * PI_)  ! 09/14/2020 BUGFIXED
      isca(1,:) = Dust_wav(iwav)             ! Wavelength (um)
      isca(2,:) = Dust_siz(:)                ! Maximum dimension
      isca(3,:) = wrks(3,:) * sfac**3.0_R_   ! Volume     (um**3)
      isca(4,:) = wrks(2,:) * sfac**2.0_R_   ! Projected area (mu**2)
      isca(5,:) = wrks(4,:) / wrks(2,:)      ! Extinction efficiency (Cext/Aproj)
      isca(6,:) = wrks(5,:) / wrks(4,:)      ! Single scattering albedo (Csca/Cext)
      isca(7,:) = wrks(6,:) / wrks(5,:)      ! asymmetry factor (g*Csca/Csca)

      ! convert phase matrix and normalization (Csca*P?? in the database)
      do ipol = 1, Dust_npol 
         pmat(1:Dust_nang,ipol,:) = wrkp(1+(ipol-1)*Dust_nang:ipol*Dust_nang,:)
      end do
      do ipol = 2, Dust_npol ! P?? --> P??/P11
         pmat(1:Dust_nang,ipol,:) = pmat(1:Dust_nang,ipol,:)/pmat(1:Dust_nang,1,:)
      end do 
      do iang = 1, Dust_nang ! Csca*P11 --> P11
         pmat(iang,1,:) = pmat(iang,1,:) / wrks(5,:) 
      end do
      !//Note (09/14/2020): TAMUdust2020 is now using a kernel technique.
      !
      ! Reference: Twomey, S. (1977). Introduction to the mathematics of
      ! inversion in remote sensing and indirect measurements. New York:
      ! Elsevier.

   end subroutine TAMUdust2020_convert


   !+
   ! TAMUdust2020: Open a database
   !-
   subroutine TAMUdust2020_database_open(iunit)

     integer,      intent(out) :: iunit(:)   ! Unit numbers
     character(len=256) :: fname(Dust_nfile)       ! Filename of input file
     integer  :: ifile, ipol


     if ( Dust_dattyp == 'DataTable' ) then

        ! Set filenames
        fname(1)  = "P11.dat"
        fname(2)  = "P12.dat"
        fname(3)  = "P22.dat"
        fname(4)  = "P33.dat"
        fname(5)  = "P43.dat"
        fname(6)  = "P44.dat"
        fname(7)  = "isca.dat"
        fname(8)  = "ReadMe.md"
        iunit(1)  = 701
        iunit(2)  = 702
        iunit(3)  = 703
        iunit(4)  = 704
        iunit(5)  = 705
        iunit(6)  = 706
        iunit(7)  = 707
        iunit(8)  = 708


        ! Open the isca.dat and P??.dat
        do ifile = 1, Dust_nfile
           call open_seq(iunit(ifile), trim(Dust_outpath) //'/'// fname(ifile), 'unknown')
        end do
        do ipol = 1, Dust_npol 
            write(iunit(ipol), '(1000es16.7e2)') Dust_ang(1:Dust_nang)
        end do

     else if ( Dust_dattyp == 'SingleCal' ) then
  
        ! Set filenames
        fname(1)  = "PMat.txt"
        fname(2)  = "isca.txt"
        iunit(1)  = 701
        iunit(2)  = 702


        ! Open the isca.dat and PMat.dat
        do ifile = 1, Dust_nfile
           call open_seq(iunit(ifile), trim(Dust_outpath) //'/'// fname(ifile), 'unknown')
        end do
        write(iunit(1), '(a112)') 'ScatteringAngle P11             P12/P11         P22/P11         P33/P11         '//&
                                  'P43/P11         P44/P11         '
        write(iunit(2), '(a112)') 'Wavelength      MaxDimension    Volume(um^3)    ProjArea(um^2)  Qext            '//&
                                  'Omega           g               '   
     end if

   end subroutine TAMUdust2020_database_open


   !+ 
   ! TAMUdust2020: Output a database 
   !-
   subroutine TAMUdust2020_write_database(iuo, isca, pmat)
 
      integer, intent(in)  :: iuo(:)      ! output file unit
      real(R_), intent(in) :: isca(:,:)   ! (nvar=7) isca.dat
      real(R_), intent(in) :: pmat(:,:,:) ! (nang,npol) pmat.dat
      integer :: isiz, ipol

      ! write 
      do isiz = 1, Dust_nsiz
         do ipol = 1, Dust_npol 
            write(iuo(ipol), '(1000es16.7e2)') pmat(:,ipol,isiz)
         end do
         write(iuo(7), '(7es16.7e2)') isca(:,isiz)
      end do

   end subroutine TAMUdust2020_write_database

   
   !+ 
   ! TAMUdust2020: Output a single case
   !-
   subroutine TAMUdust2020_write_single(iuo, isca, pmat)

      integer, intent(in)  :: iuo(:)      ! output file unit
      real(R_), intent(in) :: isca(:,:)   ! (nvar=7) isca.dat
      real(R_), intent(in) :: pmat(:,:,:) ! (nang,npol) pmat.dat
      integer :: iang, ipol

      ! write 
      do iang = 1, Dust_nang
         write(iuo(1), '(7es16.7e2)') Dust_ang(iang), (pmat(iang,ipol,1), ipol = 1, Dust_npol)
      end do
      write(iuo(2), '(7es16.7e2)') isca(:,1)

   end subroutine TAMUdust2020_write_single


   !+ 
   ! TAMUdust2020: Output a ReadMe file 
   !-
   subroutine TAMUdust2020_write_ReadMe(iuo)
 
      integer, intent(in)  :: iuo(:)      ! output file unit

      ! write 
      write(iuo(8), '(a)') &
        '------------ ReadMe.md, the ReadMe file for this database --------------'&
      //new_line('A')//new_line('A') &
      //'  1. General information'//new_line('A')//new_line('A') &                                                
      //'    This scattering property database is created from the TAMUdust2020'//new_line('A') &
      //'    database. The organization of this database is described in this'//new_line('A') &
      //'    ReadMe file.'//new_line('A')//new_line('A') &
      //'      Lead developer: Dr. Masanori Saito (masa.saito@tamu.edu)'//new_line('A') &
      //'      PI:             Dr. Ping Yang'//new_line('A')//new_line('A') &
      //'    Reference:'//new_line('A')//new_line('A') &
      //'      1). Saito, M., P. Yang, J. Ding, and X. Liu, A comprehensive'//new_line('A') &   
      //'          database of the optical properties of irregular aerosol'//new_line('A') &   
      //'          particles for radiative transfer simulations, J. Atmos.'//new_line('A') &   
      //'          Sci., submitted.'//new_line('A')//new_line('A')//new_line('A') &   
      //'  2. Database information'//new_line('A')//new_line('A') &
      //'    The database includes eight ASCII files including:'//new_line('A')//new_line('A') &
      //'      isca.dat  :: Geometric and scattering properties'//new_line('A') & 
      //'      P11.dat   :: Scattering phase fuction P11'//new_line('A') & 
      //'      P12.dat   :: Normalized scattering phase matrix element P12/P11'//new_line('A') & 
      //'      P22.dat   :: Normalized scattering phase matrix element P22/P11'//new_line('A') & 
      //'      P33.dat   :: Normalized scattering phase matrix element P33/P11'//new_line('A') & 
      //'      P43.dat   :: Normalized scattering phase matrix element P43/P11'//new_line('A') & 
      //'      P44.dat   :: Normalized scattering phase matrix element P44/P11'//new_line('A') & 
      //'      ReadMe.md :: Information of the database (THIS file).'//new_line('A')  
      write(iuo(8), '(a,i4,a,i4,a,i7,a)') &
        '    "isca.dat" consists of', &
        Dust_nwav,' x',Dust_nsiz,' =',Dust_nwav*Dust_nsiz,' lines (i.e., nwav*nsiz;'//new_line('A') &
      //'    size loops first) and 7 columns as follows:'//new_line('A')//new_line('A') & 
      //'      1). wavelength (um)' 
      if      (Dust_siztyp == 'MaximumDim') then
         write(iuo(8), '(a)') '      2). maximum dimension of a particle (um)'
      else if (Dust_siztyp == 'ProjecEqvS') then
         write(iuo(8), '(a)') '      2). projected-area-equivalent sphere diameter (um)'
      else if (Dust_siztyp == 'VolumeEqvS') then
         write(iuo(8), '(a)') '      2). volume-equivalent sphere diameter (um)'
      end if
      write(iuo(8), '(a)') &
        '      3). volume of a particle (um^3)'//new_line('A') & 
      //'      4). projected area of a particle (um^2)'//new_line('A') & 
      //'      5). extinction efficiency'//new_line('A') & 
      //'      6). single-scattering albedo'//new_line('A') & 
      //'      7). asymmetry factor'//new_line('A')//new_line('A') & 
      //'    where nwav is # of wavelength bins, and nsiz is # of size bins.'//new_line('A')
      write(iuo(8), '(a,i4,a,i4,a,i7,a)') &
        '    "P??.dat" consists of',Dust_nwav,' x',Dust_nsiz,' + 1 =',Dust_nwav*Dust_nsiz+1, &
        ' (nwav*nsiz+1) lines'//new_line('A') &
      //'    and 498 columns corresponding to the scatteringangle. The first'//new_line('A') &
      //'    line shows scattering angles (0-180 degree) and the rest of lines'//new_line('A') & 
      //'    shows P?? values for each wavelength and size (size loops first).' &
      //new_line('A')//new_line('A')//new_line('A') & 
      //'  3. Requirements'//new_line('A')//new_line('A') & 
      //'    The only requirement in regards to utilizing scattering properties'//new_line('A') &
      //'    generated from the database for research is to acknowledge our'//new_line('A') & 
      //'    contribution in a paper to be published by:'//new_line('A')//new_line('A') & 
      //'      1. citing our paper in a relevant place of main text (this will be'//new_line('A') & 
      //'         available upon the acceptance of our paper)'//new_line('A') & 
      //'      2. adding the following description in Acknowledgement section or'//new_line('A') & 
      //'         Data Availability section:'//new_line('A')//new_line('A') & 
      //'    "The scattering properties are obtained from TAMUdust2020'//new_line('A') & 
      //'    (https://sites.google.com/site/masanorisaitophd/data-and-resources)"' &
      //new_line('A')//new_line('A')//new_line('A') & 
      //'  4. Contact Information'//new_line('A')//new_line('A') & 
      //'    Dr. Masanori Saito (masa.saito@tamu.edu)'//new_line('A') &
      //'    Department of Atmospheric Sciences, Texas A&M University.'

   end subroutine TAMUdust2020_write_ReadMe


end module tamudust2020

