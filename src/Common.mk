.SUFFIXES :
.SUFFIXES : .F90 .f90 .f .o

## Compiler-specific macros:
# FC      : Fortran 90/95 compiler command
# FCFLAGS : general options for the compiler


## GNU Fortran
FC = gfortran
#FCFLAGS = -O -Wall
FCFLAGS = -Wall -fbounds-check
#LAPACKFLAGS = 
#LAPACKFLAGS = -liblapack -libblas

## G95
#FC = g95
#FCFLAGS = -O3 -Wall
#FCFLAGS = -Wall -fbounds-check

## Intel Fortran
#FC = ifort
#FCFLAGS = -O -warn
#FCFLAGS = -ipo -O
#FCFLAGS = -ipo -O3 -xP -static
# - Note: IPO option (-ipo) requires that AR = xiar. See the bottom section of this file.
#FCFLAGS = -warn -CB
#FCFLAGS = -warn -std95 -check all -CB -traceback
#-convert big_endian

## PGI Fortran
#FC = pgfortran
#FCFLAGS = -fast
#FCFLAGS = -C -traceback


### Other system-dependent macros
# - Options for files to be preprocessed (e.g., by CPP)
#CPPFLAGS = -I..
# - Archive command to make a static library (usually "ar", but use "xiar" for "ifort -ipo")
#AR = xiar
#AR = ar
# - ranlib command: Set as "" if ranlib command is not available.
#RANLIB  = ranlib
# - Option for defining a precompiler (CPP) macro
#FC_DEF = -D
# - Option for compiling MPI codes
#FC_MPI = -lmpi

# - Install path
BINDIR = ../bin

# - Library directory
LIB_DIR = ./hparx

# - Library object files
LIB_OBJS = $(LIB_DIR)/libhparx.a
# - Library module include options
LIB_FLAGS = -I$(LIB_DIR)

