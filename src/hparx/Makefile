# Makefile to make HPARX library

INC_COMMON = ./Common.mk
include $(INC_COMMON)

## Macros
TARGETS = libhparx.a
LIB_SRCS = globals.f90 hparx_base.f90 hparx_file.f90 hparx_vecmat.f90 hparx_lina.f90 \
           hparx_math.f90 hparx.f90 hparx_nons.f90
LIB_OBJS = $(LIB_SRCS:.f90=.o)
BUPIDX  = $(shell date "+%Y%m%d_%H%M")
BUPFILE = hparx-src.$(BUPIDX).tar.gz

## Commands & rules
all : $(TARGETS)

clean :
	rm -f *.o *.mod $(TARGETS) *~
echo :
	@ echo $(LIB_SRCS)
print :
	@ cat $(LIB_SRCS)
backup : clean
	@ rm -f $(BUPFILE)
	tar cvzf $(BUPFILE) *.f90 Makefile *.mk
.f90.o : 
	$(FC) -c $< $(FCFLAGS)

## Targets
libhparx.a : $(LIB_OBJS)
	@ rm -f $@
	@ echo "Archiving..." $@
	@ if test $(RANLIB); then \
		$(AR) cr $@ $^; ranlib $@; \
        else \
		$(AR) crs $@ $^; \
	fi
	@ echo "Done."

## Dependency
hparx.o : globals.o hparx_base.o hparx_vecmat.o hparx_lina.o hparx_math.o hparx_file.o hparx_nons.o
hparx_base.o   : globals.o
hparx_file.o   : globals.o hparx_base.o
hparx_vecmat.o : globals.o
hparx_lina.o   : globals.o hparx_base.o hparx_vecmat.o
hparx_math.o   : globals.o hparx_lina.o hparx_vecmat.o
hparx_nons.o   : globals.o hparx_base.o
