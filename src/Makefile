
INC_COMMON = ./Common.mk
include $(INC_COMMON)

## Macros
TARGETS  = tamudust2020create tamudust2020single
SHLS     = $(STARGETS)
BUPIDX   = $(shell date "+%Y%m%d_%H%M")
BUPFILE  = tamudust2020-src.$(BUPIDX).tar.gz

## Files
SRCS = main_TAMUdust2020_databs.f90 main_TAMUdust2020_single.f90 TAMUdust2020_module.f90
OBJS = $(SRCS:.f90=.o)

## Commands
all : $(TARGETS) $(STARGETS)

install : $(TARGETS) $(STARGETS)
	@ mkdir -p $(BINDIR)
	cp $(TARGETS) $(STARGETS) $(BINDIR)
clean  :
	rm -f *.o *.mod $(TARGETS) *~
echo   :
	@ echo $(SRCS) $(STARGETS)
print  :
	@ cat $(SRCS) $(STARGETS)
backup : $(BUPFILE)
$(BUPFILE) :
	@ rm -f $(BUPFILE)
	tar cvzf $(BUPFILE) $(SRCS) $(STARGETS) Makefile *.mk

## Suffix rules
.f90.o :
	$(FC) -c $(FCFLAGS) $(LIB_FLAGS) $(LAPACK_FLAGS) $<

## Targets
tamudust2020create   : main_TAMUdust2020_databs.o TAMUdust2020_module.o $(LIB_OBJS) 
	$(FC) -o $@ $(FCFLAGS) $^
tamudust2020single   : main_TAMUdust2020_single.o TAMUdust2020_module.o $(LIB_OBJS) 
	$(FC) -o $@ $(FCFLAGS) $^

## Dependency
main_TAMUdust2020_databs.o        : $(LIB_OBJS) TAMUdust2020_module.o
main_TAMUdust2020_single.o        : $(LIB_OBJS) TAMUdust2020_module.o
TAMUdust2020_module.o             : $(LIB_OBJS) 
