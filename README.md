<<<<<<< HEAD
------------- READ-ME file for the TAMUdust2020 database ---------------
=======
------READ-ME file for the TAMUdust2020 database------
>>>>>>> c1312797ff9e90923ccd835d56677b1619403258

 This is a flexible user-friendly database for optical properties of 
 dust aerosol and volcanic ash particles. Various ensembles of 20 
 irregular hexahedral particle models are applied in this database. 
 The organization of this database is described in this ReadMe file.


 Lead developer: Dr. Masanori Saito (masa.saito@tamu.edu)
 PI: 		 Dr. Ping Yang 


 History
 (MM/DD/YYYY)	(Author)	(Contents)
 06/01/2018	Masa Saito	Project initiated.
 03/01/2019 	Masa Saito	Extended the applicability from 3 lidar
 				wavelengths to shortwave (0.2 to 4 
                                micron).
 03/16/2020	Masa Saito	Sphericity variation is considered.
 05/01/2020	Masa Saito	Extended the applicability from 
                                shortwave to the entire wavelengths 
                                (0.2-200 micron).
 07/30/2020	Masa Saito	The database Version 0.9 has been 
                                completed. 
<<<<<<< HEAD
 10/07/2020     Masa Saito	Source code published.
 04/22/2021     Masa Saito	Data kernel published. 
                                Version 1.0 released.
=======
 09/30/2020     Masa Saito      The database Version 0.9.9 has been 
                                completed (This will automatically 
                                be Version 1.0 after publication). 
 10/07/2020     Masa Saito      Source code published.
  
 
 ***Note: Version 1.0 will be available after the publication of the 
          paper on the TAMUdust2020 database. 
>>>>>>> c1312797ff9e90923ccd835d56677b1619403258



1. General information

   The TAMUdust2020 database is designed for users who would like to 
   obtain the single-scattering properties of a user-defined aerosol 
   optical property model. For example, if users have the data of 
   wavelength-dependent complex refractive of aerosol particles, they 
   can create the single-scattering properties of the particles. 

   Reference:

     1). Saito, M., P. Yang, J. Ding, and X. Liu, A comprehensive 
         database of the optical properties of irregular aerosol 
         particles for radiative transfer simulations, J. Atmos. Sci., 
<<<<<<< HEAD
         in press, https://doi.org/10.1175/JAS-D-20-0338.1.
=======
         submitted.
>>>>>>> c1312797ff9e90923ccd835d56677b1619403258



   1.1 Code and Database

     The TAMUdust2020 database consists of the following directories:
     
       src		:: Source codes
       src/hparx  	:: HPARX library used for the TAMUdust2020 
                           database creation
       bin		:: Executable command path
       examples		:: Example files
<<<<<<< HEAD
       examples/output  :: Directory for the output from the test run 
=======
>>>>>>> c1312797ff9e90923ccd835d56677b1619403258
       examples/params  :: Directory for the parameter files for the 
                           test run 

     In addition to these three directories, users will need to obtain 
<<<<<<< HEAD
     the following directory data from Zenodo.org (see the database 
     webpage: https://sites.google.com/site/masanorisaitophd/
     data-and-resources/tamudust2020).
=======
     the following directory data from:

     (Table path will be publicly available after the publication of 
      the paper. Please contact to the lead developer to obtain the 
      path for the database.)
>>>>>>> c1312797ff9e90923ccd835d56677b1619403258

       tables	 	:: Data kernel
       tables/SW 	:: Data kernel for shortwave (0.2-4 um)
       tables/LW 	:: Data kernel for longwave (>4 um)
<<<<<<< HEAD
      
     The latest version of the source codes and examples may be 
     available from https://github.com/masasaito/TAMUdust2020.
=======
>>>>>>> c1312797ff9e90923ccd835d56677b1619403258



   1.2 Variable Definition

      W   :: wavelength (micron, um)
      D   :: maximum dimension (um), defined as an average of
             circumscribed diameters of individual particles.
      V   :: volume of a particle (um^3)  
      A   :: projected area (um^2)
      Qe  :: extinction efficiency
      SSA :: single-scattering albedo
      g   :: asymmetry factor
      x   :: size parameters defined as x = pi*D/W
      mr  :: real part of refractive index
      mi  :: imaginary part of refractive index
      S   :: sphericity
      P?? :: phase matrix elements (?? includes 11, 12, 22, 33, 43, 
             and 44)

      Whole steradian integral of P11 is 4*pi:
    
           pi     2*pi
        int    int     P11 d(zenith angle)d(azimuth angle)=4*pi
           0      0

      Other independent phase matrix elements P12, P22, ..., and P44 
      are normalized by P11.


      The database includes two kernels (SW and LW), and each kernel 
      covers ranges of optical/physical properties differently as:

      SW Kernel			LW Kernel           
      x   :: <=11800		x   :: <=1470
      mr  :: 1.37--1.7		mr  :: 0.4--3.2
      mi  :: 0.0001--0.1 	mi  :: 0.001--4.0
      S   :: 0.695--0.785	S   :: 0.695--0.785
         


 2. Installation

   First, to install the database creation codes, edit the following 
   files:

     ./src/Common.mk
     ./src/hparx/Common.mk  
     
   In these files, users can specify their compiler (FC; e.g., gfortran)
   and associate options (FCFLAGS).

   Then, command as follows:

     $ cd ./src/hparx
     $ make
     $ cd ../
     $ make install
 
   In the end, the following two executable files will be created in the 
   "bin" directory:

     tamudust2020single	:: Create the single-scattering properties 
     			   of a user-defined aerosol particle.
     tamudust2020create :: Create the single-scattering property 
  			   database of user-defined aerosol optical
  			   property model.



      
 3. How to run the codes

   Go to the example directory (examples), and users will find several
   files and directories as follows:
 
     job_runTEST.sh		 :: job script for test run
<<<<<<< HEAD
     output/ 			 :: a directory for output 
=======
>>>>>>> c1312797ff9e90923ccd835d56677b1619403258
     params/			 :: a directory containing parameter 
     				    files
     TAMUdust2020create_exp1.nml :: namelist for database creation for 
				    SW domain 
     TAMUdust2020create_exp2.nml :: namelist for database creation for 
				    LW domain
     TAMUdust2020single_exp1.nml :: namelist for single-scattering 
  				    property creation for SW domain
     TAMUdust2020single_exp2.nml :: namelist for single-scattering 
				    property creation for LW domain

   Simply, type the following command

     $ ./job_runTEST.sh
 
   Then, users will notice that the output files are created at
<<<<<<< HEAD
   ./output/{DIRECTORY_FOR_EACH_EXAMPLE}/.  
=======
   the 'example' directory.  
>>>>>>> c1312797ff9e90923ccd835d56677b1619403258



   3.1 How to run "tamudust2020single"

     This code creates the single-scattering properties of an 
     aerosol model defined by a user with the following parameters:

      D  :: maximum diameter (um)
      W  :: wavelengths (um)
      mr :: real part of refractive index
      mi :: imaginary part of refractive index
      S  :: sphericity

     The namelist for tamudust2020single run contains (e.g., exp1):

        ! path
        TAMUdust2020path = '../tables/SW'           !(FIXED)  
        angfile    = './params/TAMUdust2020_Angle'  !(FIXED)

        ! table configuration
        nsiztab   = 169                             !(FIXED)
        nrefretab = 8                               !(FIXED)       
        nrefimtab = 31                              !(FIXED)
        nsphtab   = 6                               !(FIXED)

     Users do not have to change the table configurations that are fixed
     for the corresponding data kernel and angle file that is fixed.

     To run tamudust2020single, simply type:

       $ ../bin/tamudust2020single {NAMELIST} D W mr mi S {OUTPUT_PATH}



   3.2 How to run "tamudust2020create"

     This code creates the single-scattering property database of an 
     aerosol model specified in a namelist (e.g., exp1):

       ! path
       TAMUdust2020path = '../tables/SW'                         !(FIXED)
<<<<<<< HEAD
       outpath    = './output/TAMUdust2020create_exp1'           !(USER)   
=======
       outpath    = './output_folder'                            !(USER)   
>>>>>>> c1312797ff9e90923ccd835d56677b1619403258
       angfile    = './params/TAMUdust2020_Angle'                !(FIXED)
       wavfile    = './params/TAMUdust2020_Wavelength_exp1'      !(USER)
       sizfile    = './params/TAMUdust2020_Size'                 !(USER)
       refidxfile = './params/TAMUdust2020_RefractiveIndex_exp1' !(USER)

       ! table configuration                                     
       nsiztab   = 169                                           !(FIXED)
       nrefretab = 8                                             !(FIXED)
       nrefimtab = 31                                            !(FIXED)
       nsphtab   = 6                                             !(FIXED)

       ! database configuration
       nsiz = 43              ! # of maximum dimension bin       !(USER)
       nwav = 39              ! # of wavelength bin              !(USER) 
       spher = 7.0000000E-01  ! Sphericity of particles          !(USER)

     Users do not have to change the table configurations that are fixed
     for the corresponding data kernel and angle file that is fixed. 

     For other files/parameters marked as (USER) is a user-defined one 
<<<<<<< HEAD
     editable as explained in the subsections below, including 
     wavelengths, size, refractive index, and sphericity.
=======
     editable as explained in the subsections below, including wavelengths, 
     size, refractive index, and sphericity.
>>>>>>> c1312797ff9e90923ccd835d56677b1619403258

     To run tamudust2020create, simply type:

       $ ../bin/tamudust2020create {NAMELIST}


     3.2.1 Wavelength specification

       'wavfile' contains wavelengths in lines 1 to 'nwav'. For example,
       in the case of nwav=7, wavfile contains:

       2.0000000e-01
       3.0000000e-01
       4.0000000e-01
       5.0000000e-01
       6.0000000e-01
       7.0000000e-01
       8.0000000e-01
       

     3.2.2 Particle size specification

       'sizfile' contains particle maximum dimension in lines 1 to 
       'nsiz'. For example, in the case of nsiz=9, sizfile contains:

       1.0000000e-02
       3.2000000e-02
       1.0000000e-01
       3.2000000e-01
       1.0000000e+00
       3.2000000e+00
       1.0000000e+01
       3.2000000e+01
       1.0000000e+02


     3.2.3 Refractive index specification

       'refidxfile' contains spectral refractive index in lines 1 to 
       'nwav'. The first and second columns correspond to the real and
       imaginary part of the refractive index. For example, in the case 
       of nwav=7, refidxfile contains:

       1.56751E+00      1.43611E-01
       1.55622E+00      6.45710E-02
       1.55967E+00      3.94937E-02
       1.55236E+00      3.34802E-02
       1.54350E+00      2.86088E-02
       1.52884E+00      3.01948E-02
       1.51796E+00      3.07916E-02


     3.2.4 Sphericity specification
      
       'spher' indicate user-defined sphericity. 



 
4. Output 

   4.1 "tamudust2020single"

     Output is two ASCII files including:

       isca.txt	:: Geometric and scattering properties 
       PMat.txt :: Scattering phase matrix elements

    'isca.txt' consists of a single line and 7 columns as follows:

    1). wavelength (um) 
    2). maximum dimension of the particle (um)
    3). volume of particle (um^3) 
    4). projected area (um^2)
    5). extinction efficiency
    6). single-scattering albedo 
    7). asymmetry factor

    'PMat.txt' consists of 498 lines corresponding to the scattering 
    angle and 7 columns as follows:

    1). scattering angle (degree)
    2). P11
    3). P12/P11
    4). P22/P11
    5). P33/P11
    6). P43/P11
    7). P44/P11



   4.2 "tamudust2020create"

     Output is seven ASCII files including:

       isca.dat	:: Geometric and scattering properties 
       P11.dat  :: Scattering phase matrix element P11
       P12.dat  :: Normalized scattering phase matrix element P12/P11
       P22.dat  :: Normalized scattering phase matrix element P22/P11
       P33.dat  :: Normalized scattering phase matrix element P33/P11
       P43.dat  :: Normalized scattering phase matrix element P43/P11
       P44.dat  :: Normalized scattering phase matrix element P44/P11.

    'isca.dat' consists of nwav*nsiz lines (size loops first) and 
    7 columns as follows:
     
    1). wavelength (um) 
    2). maximum dimension of the particle (um)
    3). volume of particle (um^3) 
    4). projected area (um^2)
    5). extinction efficiency
    6). single-scattering albedo 
    7). asymmetry factor

    For example, 'isca.dat' with nwav=10 and nsiz=30 consists of 
    10*30=300 lines in which the first 30 lines are for different 
    sizes with a first wavelength.


    'P??.dat' consists of nwav*nsiz+1 lines and 498 columns corresponding
    to the scattering angle. The first line shows scattering angles 
    (0-180 degree) and the rest of lines shows P?? values for each 
    wavelength and size (size loops first).


  

  5. Requirements

    The lead developer, Dr. Masanori Saito encourages the research
    community to utilize the TAMUdust2020 database for aerosol studies. 
    The only requirement in regards to utilizing scattering properties 
    generated from the database for research is to acknowledge our 
    contribution in a paper to be published by:

<<<<<<< HEAD
      1. citing Saito et al. (2021) in a relevant section of the main 
         text. 
      2. adding the following description in Acknowledgement section 
         or Data Availability section:
         "The scattering properties are obtained from TAMUdust2020."
=======
      1. citing our paper in a relevant place of main text (this will be 
         available upon the acceptance of our paper)
      2. adding the following description in Acknowledgement section or 
         Data Availability section:
      
    "The scattering properties are obtained from TAMUdust2020
    (https://github.com/masasaito/TAMUdust2020)."

   ***Note: This will be modified after official publication of our paper.
>>>>>>> c1312797ff9e90923ccd835d56677b1619403258
      



  6. Contact Information

     Dr. Masanori Saito (masa.saito@tamu.edu) 
     Department of Atmospheric Sciences, Texas A&M University. 

 
