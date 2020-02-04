SHELL=	/bin/sh

include $(ESMFMKFILE)

CMD= chgres_cube.exe

OBJS = chgres.o \
       atmosphere.o \
       input_data.o \
       surface.o \
       model_grid.o \
       program_setup.o \
       write_data.o \
       search_util.o \
       static_data.o \
       utils.o \
       grib2_util.o

$(CMD):	$(OBJS)
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -o $(CMD) $(OBJS) $(SP_LIBd) $(NEMSIO_LIB) $(BACIO_LIB4) $(W3NCO_LIBd) $(SFCIO_LIB4) $(SIGIO_LIB4) $(ESMF_F90LINKPATHS) $(ESMF_F90ESMFLINKRPATHS) $(ESMF_F90ESMFLINKLIBS) $(WGRIB2_LIB)

model_grid.o:  program_setup.o model_grid.F90
	$(FCOMP) $(FFLAGS) -I$(NEMSIO_INC) -I$(SFCIO_INC4) -I$(SIGIO_INC4) -I$(WGRIB2API_INC) $(ESMF_F90COMPILEPATHS) -c model_grid.F90

utils.o:  utils.f90
	$(FCOMP) $(FFLAGS) -c $(ESMF_F90COMPILEPATHS) utils.f90

grib2_util.o:  model_grid.o program_setup.o grib2_util.F90
	$(FCOMP) $(FFLAGS) -c $(ESMF_F90COMPILEPATHS) grib2_util.F90

program_setup.o:  program_setup.f90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c program_setup.f90

atmosphere.o:  program_setup.o model_grid.o input_data.o atmosphere.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c atmosphere.F90

chgres.o: atmosphere.o model_grid.o program_setup.o surface.o chgres.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c chgres.F90

write_data.o: atmosphere.o  model_grid.o program_setup.o surface.o static_data.o write_data.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c write_data.F90

input_data.o: program_setup.o model_grid.o grib2_util.o input_data.F90
	$(FCOMP) $(FFLAGS) -I$(NEMSIO_INC) -I$(SIGIO_INC4) -I$(SFCIO_INC4) -I$(WGRIB2API_INC) $(ESMF_F90COMPILEPATHS) -c input_data.F90

surface.o: search_util.o model_grid.o input_data.o program_setup.o static_data.o surface.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c surface.F90

search_util.o: search_util.f90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c search_util.f90

static_data.o: model_grid.o program_setup.o static_data.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c static_data.F90

clean:
	rm -f *.o *.mod ${CMD} ../../exec/${CMD}

install:
	mv $(CMD) ../../exec
