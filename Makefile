COMPILE.f90 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -c
MAKEMOD.f90 = $(FC) $(FCFLAGS) $(TARGET_ARCH) -fsyntax-only -c

MODS = $(wildcard module*.f90)
MOD_OBJS = $(patsubst %.f90,%.o,$(MODS))

MAIN_SOURCE = main.f90
MAIN_OBJ = $(patsubst %.f90,%.o,$(MAIN_SOURCE))

FC=gfortran
FLFLAGS = -g -fbounds-check -fbacktrace -fcheck=all -ffree-line-length-0
FFLAGSS = -Wall
PROGRAM=CCL1D

.SUFFIXES: .o .f90 .cc

.PHONY: debug default clean

.cc.o:
	$(FC) $(FLFLAGS) -I $(CDIR) -c $<

.f90.o:
	$(FC) -c $(FLFLAGS) $<

CCL1D: $(MOD_OBJS) $(MAIN_OBJ)
	$(FC) $(FFLAGSS) -o CCL1D $(MOD_OBJS) $(MAIN_OBJ)

clean:
	rm -f CCL1D *.o *.mod *.f90~ core

#Dependencies
CCL1D.o: module_data.o module_imexport.o module_geometry.o module_solvers.o module_eos.o
module_geometry.o: module_data.o
module_imexport.o: module_data.o module_geometry.o module_eos.o
module_solvers.o: module_data.o module_eos.o
