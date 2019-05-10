SHELL=	/bin/sh

include $(ESMFMKFILE)

CMD= sfc_climo_gen

OBJS = driver.o \
       model_grid.o \
       output.o \
       program_setup.o \
       source_grid.o \
       utils.o \
       search.o \
       interp.o

$(CMD):	$(OBJS)
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -o $(CMD) $(OBJS) $(ESMF_F90LINKPATHS) $(ESMF_F90LINKRPATHS) $(ESMF_F90ESMFLINKLIBS) -g -traceback

interp.o:  source_grid.o model_grid.o utils.o interp.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c interp.F90

utils.o:  utils.f90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c utils.f90

program_setup.o:  program_setup.f90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c program_setup.f90

model_grid.o:  utils.o program_setup.o model_grid.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c model_grid.F90

source_grid.o:  utils.o source_grid.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c source_grid.F90

search.o:  search.f90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c search.f90

output.o: program_setup.o model_grid.o source_grid.o utils.o  output.f90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c output.f90

driver.o: model_grid.o source_grid.o program_setup.o utils.o driver.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c driver.F90

install:
	mv $(CMD) ../../exec
clean:
	rm -f *.o *.mod ${CMD} ../../exec/${CMD}
