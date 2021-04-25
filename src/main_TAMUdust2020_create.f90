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
! 
!-
program tamudust2020create

  use globals
  use hparx, only : open_dir, gridIdx_loc, close_all, open_seq, getCmdArgs, err_read, intpTab_Akima_1R
  implicit none

  ! TAMUdust2020 binary kernel table
  integer, save      :: Tab_nsiz   = 169      ! # of size bin
  integer, parameter :: Tab_nang   = 498      ! # of scattering angle bin
  integer, parameter :: Tab_npol   = 6        ! # of polarization elements
  integer, parameter :: Tab_nmet   = 9        ! # of meta information
  integer, save      :: Tab_nsph   = 6        ! # of sphericity bin
  integer, save      :: Tab_nReIdx = 8        ! # of real part of dust refractive index bin
  integer, save      :: Tab_nImIdx = 31       ! # of imaginary part of dust reflactive index bin
  integer, save      :: Tab_nref              ! # of refractive index bin
  integer, save      :: Tab_ReIdxSoft_ID = 1  ! ID for optically soft refractive index (=1.1)
  integer, parameter :: Tab_nfile  = 2        ! binary pmat and isca files
  real(R4_), allocatable :: TAMUdust2020STab(:,:,:,:) ! isca.dat kernel table
  real(R4_), allocatable :: TAMUdust2020PTab(:,:,:,:) ! Phase matrix dat kernel table 
  real(R_), allocatable  :: Tab_siz(:), Tab_ReIdx(:), Tab_ImIdx(:), Tab_sph(:) 

  ! TAMUdust2020 ASCII user-defined aerosol scattering property database
  character(len=256)  :: Dust_outpath       ! output dust database path
  integer, parameter  :: Dust_nfile = 7     ! # of files
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

  ! TAMUdust2020 Core Process
  call TAMUdust2020__init()
  call TAMUdust2020__create()
  call TAMUdust2020__final()

stop
contains


   !+
   ! TAMUdust2020: Read namelist
   !-
   subroutine TAMUdust2020__init()
     
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
      character(len=256)       :: TAMUdust2020path  ! TAMUdust2020 binary table path
      integer,  save           :: nsiz, nwav        ! # of size/wavelength bins
      integer, save            :: nsiztab, nrefretab, nrefimtab, nsphtab
      real(R_), save           :: spher
      namelist /TAMUdust2020/TAMUdust2020path, outpath, sizfile, wavfile, refidxfile, angfile, &
                             nsiz, nwav, spher, nsiztab, nrefretab, nrefimtab, nsphtab            

      argmsg = ' Usage: tamudust2020create nmlst > output '  
      narg = 1
      call getCmdArgs(narg, argv, argmsg)
      open(IUI, file = trim(adjustl(argv(1))), status = 'old')
      read(IUI, nml=TAMUdust2020, iostat=ios)
      call err_read(ios, IUI, 'Main: Error in reading the namelist data.')
      close(IUI)

 
      ! Initialization
      Dust_outpath = outpath
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


      ! Get meta information
      Tab_ReIdxSoft_ID = 1
      Tab_siz(1:Tab_nsiz) = TAMUdust2020STab(1,1,1:Tab_nsiz,1)
      do ire = 1, Tab_nReIdx
         Tab_ReIdx(ire) =  TAMUdust2020STab(7,1,1,1+(ire-1)*Tab_nImIdx)
         if (Tab_ReIdx(ire) < 1.1005_R_) Tab_ReIdxSoft_ID = ire
         do iim = 1, Tab_nImIdx
            Tab_ImIdx(iim) =  TAMUdust2020STab(8,1,1,iim)
         end do
      end do
      do isph = 1, Tab_nsph
         Tab_sph(isph) = TAMUdust2020STab(9,isph,1,1)
      end do 


      ! Warnings
      if (maxval(Dust_siz)*2.0_R_*PI_/minval(Dust_wav) > maxval(Tab_siz)) &
      print*, 'Warning: large sizes in the database may be based on extrapolation. Size parameter ', &
      maxval(Dust_siz)*2.0_R_*PI_/minval(Dust_wav), ' is beyond the prescribed range.' 
      if (minval(Dust_refidx(:,1)) < minval(Tab_ReIdx(:)) .or. maxval(Dust_refidx(:,1)) > maxval(Tab_ReIdx(:))) &
      print*, 'Warning: Specified real part of refractive index is beyond the prescribed range. '
      if (minval(Dust_refidx(:,2)) < minval(Tab_ImIdx(:))  .or. maxval(Dust_refidx(:,2)) > maxval(Tab_ImIdx(:))) &
      print*, 'Warning: Specified imaginary part of refractive index is beyond the prescribed range. '
      if (Dust_sph < minval(Tab_sph(:))  .or. Dust_sph > maxval(Tab_sph(:))) &
      print*, 'Warning: Specified sphericity ', Dust_sph, ' is beyond the prescribed range. '


   end subroutine TAMUdust2020__init


   !+
   ! TAMUdust2020: Finalize the database
   !-
   subroutine TAMUdust2020__final()

      deallocate( Dust_siz, Dust_refidx, Dust_ang, Dust_wav )
      deallocate( TAMUdust2020STab, TAMUdust2020PTab, Tab_siz, Tab_ReIdx, Tab_ImIdx, Tab_sph ) 

   end subroutine TAMUdust2020__final


   !+
   ! TAMUdust2020: Create a database based on user-defined configuration 
   !-
   subroutine TAMUdust2020__create()
     
      integer  :: iwav, isiz, imet, ire, iim, isp
      real(R_) :: sre, sim, ssp 
      integer,  allocatable :: iuo(:)        ! I/O unit
      real(R_), allocatable :: wrk1s(:,:), wrk1p(:,:)      ! wrk for isca/p??.dat.    
      real(R_), allocatable :: wrk2s(:,:), wrk2p(:,:)      ! wrk for isca/p??.dat.    
      real(R_), allocatable :: isca(:,:),  pmat(:,:,:)     ! isca/p??.dat.    
      real(R_), allocatable :: spara(:)

      ! Set properties
      allocate( iuo(Dust_nfile), spara(Dust_nsiz) )
      allocate( wrk1s(Tab_nmet,Tab_nsiz),  wrk1p(Tab_nang*Tab_npol,Tab_nsiz)  )
      allocate( wrk2s(Tab_nmet,Dust_nsiz), wrk2p(Tab_nang*Tab_npol,Dust_nsiz) )
      allocate( isca(Dust_nmet,Dust_nsiz), pmat(Dust_nang, Dust_npol, Dust_nsiz) )
      call TAMUdust2020_database_open(iuo)

      ! Interpolation
      do iwav = 1, Dust_nwav

         ! Pickup data based on a specified refractive index
         if      (Dust_refidx(iwav,1) > 1.1005_R_) then
            call gridIdx_loc(log(Tab_ReIdx(Tab_ReIdxSoft_ID:)-1.0_R_), log(Dust_refidx(iwav,1)-1.0_R_), ire, sre, 1)
            ire = ire + Tab_ReIdxSoft_ID - 1
         else if (Dust_refidx(iwav,1) < 0.8995_R_) then
            call gridIdx_loc(log(1.0_R_-Tab_ReIdx(:)), log(1.0_R_-Dust_refidx(iwav,1)), ire, sre, 1)
            sre = 1.0_R_ - sre ! 
         else
            call gridIdx_loc(Tab_ReIdx(:), Dust_refidx(iwav,1), ire, sre, 1)
         end if
         call gridIdx_loc(log(Tab_ImIdx(:)), log(Dust_refidx(iwav,2)), iim, sim, 1)
         call gridIdx_loc(log(Tab_sph(:)),   log(Dust_sph),            isp, ssp, 1)
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
         if (minval(spara(:)) >= minval(Tab_siz(:))) then
            isiz = 1
         else

            ! Single-scattering properties: Rayleigh scattering approximation
            loop_minsiz: do isiz = 1, Tab_nsiz  
               if (spara(isiz) >= minval(Tab_siz(:))) exit loop_minsiz
            end do loop_minsiz          
            call TAMUdust2020_RayleighScat(spara(1:isiz-1), wrk1s(:,1), wrk2s(:,1:isiz-1), wrk2p(:,1:isiz-1))
         end if


         ! Single-scattering properties: Akima interpolation from the table
         do imet = 1, Tab_nmet
            wrk2s(imet,isiz:) = intpTab_Akima_1R(log(Tab_siz(:)), wrk1s(imet,:), log(spara(isiz:)), 1)  
         end do    
         do imet = 1, Tab_nang*Tab_npol
            wrk2p(imet,isiz:) = intpTab_Akima_1R(log(Tab_siz(:)), wrk1p(imet,:), log(spara(isiz:)), 1)  
         end do    

         ! write 
         call TAMUdust2020_convert(wrk2s, wrk2p, iwav, isca, pmat)
         call TAMUdust2020_write_database(iuo, isca, pmat)
      end do

      call close_all(iuo)
      deallocate( iuo, spara, wrk1s, wrk1p, wrk2s, wrk2p, isca, pmat)

   end subroutine TAMUdust2020__create


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
        sfac = spara(isiz)/Tab_siz(1)
        mref = dcmplx(wrk1s(7), -wrk1s(8))
        mre  = real(((mref**2-1.0_R_)/(mref**2+2.0_R_))**2)
        mim  = -aimag((mref**2-1.0_R_)/(mref**2+2.0_R_))
        rva  = 0.75_R_ * sfac * wrk1s(3) / wrk1s(2)              ! V/A-equivalent sphere radius (micron)
        qsca = 2.66666666_R_ * rva**4 * mre
        qabs = 4.0_R_ * rva * mim

        ! isca.dat
        wrk2s(1,isiz)   = wrk1s(1)*sfac                      ! Dmax   (Dimensionless)
        wrk2s(2,isiz)   = wrk1s(2)*sfac**2                   ! Aproj  (Dimensionless)
        wrk2s(3,isiz)   = wrk1s(3)*sfac**3                   ! Vol    (Dimensionless)
        wrk2s(4,isiz)   = (qsca + qabs)*wrk1s(2)*sfac**2     ! Cext   (Dimensionless)
        wrk2s(5,isiz)   = qsca*wrk1s(2)*sfac**2              ! Csca   (Dimensionless)
        wrk2s(6,isiz)   = 0.0_R_                             ! Csca*g (Dimensionless)
        wrk2s(7:9,isiz) = wrk1s(7:9)                         ! RefRe, RefIm, Sphericity

        ! P??.dat
        wrk2p(           1:  Tab_nang,isiz) =  0.75_R_*(1.0_R_+dcos(Dust_ang(1:Tab_nang)*DTOR_)**2) ! P11
        wrk2p(  Tab_nang+1:2*Tab_nang,isiz) = -0.75_R_*(1.0_R_-dcos(Dust_ang(1:Tab_nang)*DTOR_)**2) ! P12
        wrk2p(2*Tab_nang+1:3*Tab_nang,isiz) =  wrk2p(           1:  Tab_nang,isiz)                  ! P22
        wrk2p(3*Tab_nang+1:4*Tab_nang,isiz) =  1.5_R_ * dcos(Dust_ang(1:Tab_nang)*DTOR_)            ! P33
        wrk2p(4*Tab_nang+1:5*Tab_nang,isiz) =  0.0_R_                                               ! P43
        wrk2p(5*Tab_nang+1:6*Tab_nang,isiz) =  wrk2p(3*Tab_nang+1:4*Tab_nang,isiz)                  ! P44
        wrk2p(           1:6*Tab_nang,isiz) =  wrk2p(           1:6*Tab_nang,isiz) * qsca*wrk1s(2)*sfac**2 ! P?? --> Csca*P??
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
      !//NOTE wrks: first column (spara, Aproj, Vol, Cext, Csca, g*Csca, RefRe, RefIm)
      ! for Aproj and Vol, scaling from size parameter to actual size is needed.

      ! convert optical properties
      !sfac      = Dust_wav(iwav) / PI_
      sfac      = Dust_wav(iwav) / (2.0_R_ * PI_)  ! 09/14/2020 BUGFIXED
      isca(1,:) = Dust_wav(iwav)        ! Wavelength (um)
      isca(2,:) = Dust_siz(:)           ! Maximum dimension
      isca(3,:) = wrks(3,:) * sfac**3   ! Volume     (um**3)
      isca(4,:) = wrks(2,:) * sfac**2   ! Projected area (mu**2)
      isca(5,:) = wrks(4,:) / wrks(2,:) ! Extinction efficiency (Cext/Aproj)
      isca(6,:) = wrks(5,:) / wrks(4,:) ! Single scattering albedo (Csca/Cext)
      isca(7,:) = wrks(6,:) / wrks(5,:) ! asymmetry factor (g*Csca/Csca)

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
   ! TAMUdust2020: Open a database
   !-
   subroutine TAMUdust2020_database_open(iunit)

     integer,      intent(out) :: iunit(:)   ! Unit numbers
     character(len=256) :: fname(Dust_nfile)       ! Filename of input file
     integer  :: ifile, ipol

     ! Set filenames
     fname(1)  = "P11.dat"
     fname(2)  = "P12.dat"
     fname(3)  = "P22.dat"
     fname(4)  = "P33.dat"
     fname(5)  = "P43.dat"
     fname(6)  = "P44.dat"
     fname(7)  = "isca.dat"
     iunit(1)  = 701
     iunit(2)  = 702
     iunit(3)  = 703
     iunit(4)  = 704
     iunit(5)  = 705
     iunit(6)  = 706
     iunit(7)  = 707


     ! Open the isca.dat adn P??.dat
     do ifile = 1, Dust_nfile
        call open_seq(iunit(ifile), trim(Dust_outpath) //'/'// fname(ifile), 'unknown')
     end do
     do ipol = 1, Dust_npol 
         write(iunit(ipol), '(1000es16.7e2)') Dust_ang(1:Dust_nang)
     end do

   end subroutine TAMUdust2020_database_open


end program tamudust2020create

