.SUFFIXES :
.SUFFIXES : .F90 .f90 .f .o

### Compiler-specific macros:
# FC      : Fortran 90/95 compiler command
# FCFLAGS : general options for the compiler

## GNU Fortran
FC = gfortran
FCFLAGS = -O -Wall
#FCFLAGS = -Wall -fbounds-check
#FCFLAGS = -Wall -pedantic -std=f95 -Wuninitialized -fbacktrace -fbounds-check

## G95
#FC = g95
#FCFLAGS = -O -Wall
#FCFLAGS = -Wall -fbounds-check

## Intel Fortran
#FC = ifort
#FCFLAGS = -O -warn
#FCFLAGS = -ipo -O
#FCFLAGS = -ipo -O3 -xP -static
# - Note: IPO option (-ipo) requires that AR = xiar. See the bottom section of this file.
#FCFLAGS = -warn
#FCFLAGS = -warn -CB
#FCFLAGS = -warn -std95 -check all
#FCFLAGS = -warn -std95 -check all -CB -traceback
# - Note: Try the -heap-arrays option if you get "Segmentation fault".
#-heap-arrays
#-convert big_endian

## PGI Fortran
#FC = pgfortran
#FCFLAGS = -O
#FCFLAGS = -C -traceback


### Other system-dependent macros
# - Linker options
#LDFLAGS = -L/usr/local/lib/
# - Options for files to be preprocessed (e.g., by CPP)
#CPPFLAGS = -I..
# - Libraries to be linked
#LIBS = -lm
# - Archive command to make a static library (usually "ar", but use "xiar" for "ifort -ipo")
#AR = xiar
AR = ar
# - ranlib command: Set as "" if ranlib command is not available.
RANLIB  = ranlib
# - Option for defining a precompiler (CPP) macro
FC_DEF = -D
# - Option for compiling MPI codes
FC_MPI = -lmpi
# - Library directory
#LIB_DIR = ../hparx
# - Library object files
#LIB_OBJS = $(LIB_DIR)/libhparx.a
# - Option for adding a directory with module files (*.mod)
#FC_MOD = -I$(LIB_DIR)
