#!/bin/bash

# TAMUdust2020create_exp1 (an example dust aerosol single-cattering property database in wavelengths 0.2 to 4   µm)
../bin/tamudust2020create TAMUdust2020create_exp1.nml 

# TAMUdust2020create_exp2 (an example dust aerosol single-scattering property database in wavelengths 4.1 to 100 µm)
../bin/tamudust2020create TAMUdust2020create_exp2.nml

# TAMUdust2020single_exp1 (single-scatteing property of a given dust properties in shortwave 0.2-4 µm)  
<<<<<<< HEAD
../bin/tamudust2020single TAMUdust2020single_exp1.nml 5.0 0.532 1.5 0.003 0.7 ./output/TAMUdust2020single_exp1

# TAMUdust2020single_exp2 (single-scatteing property of a given dust properties in shortwave 4-100 µm)  
../bin/tamudust2020single TAMUdust2020single_exp2.nml 5.0 11.0  1.5 0.3   0.7 ./output/TAMUdust2020single_exp2
=======
../bin/tamudust2020single TAMUdust2020single_exp1.nml 5.0 0.532 1.5 0.003 0.7 ./

# TAMUdust2020single_exp2 (single-scatteing property of a given dust properties in shortwave 4-100 µm)  
../bin/tamudust2020single TAMUdust2020single_exp2.nml 5.0 11.0  1.5 0.3   0.7 ./
>>>>>>> c1312797ff9e90923ccd835d56677b1619403258


