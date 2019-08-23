SHELL = /bin/sh
#
# Include file specified by the environment variable ESMFMKFILE.
#
$(info )
$(info Including in the makefile the file specified by the variable ESMFMKFILE, where:)
$(info )
$(info ESMFMKFILE = "$(ESMFMKFILE)")
$(info )
include $(ESMFMKFILE)
#
# Print out variables defined in the file ESMFMKFILE (included above) 
# and used in the link step below.
#
$(info )
$(info Values of variables in file ESMFMKFILE that are used in the link step in this makefile are:)
$(info )
$(info ESMF_F90COMPILEPATHS = "$(ESMF_F90COMPILEPATHS)")
$(info ESMF_F90LINKPATHS = "$(ESMF_F90LINKPATHS)")
$(info ESMF_F90ESMFLINKPATHS = "$(ESMF_F90ESMFLINKPATHS)")
$(info ESMF_F90ESMFLINKLIBS = "$(ESMF_F90ESMFLINKLIBS)")
$(info WGRIB2_LIB = "$(WGRIB2_LIB)")
$(info WGRIB2API_LIB = "$(WGRIB2API_LIB)")
$(info )
#
# Full paths to static library files. 
#
LIBS = $(NEMSIO_LIB) \
       $(W3NCO_LIBd) \
       $(BACIO_LIB4) \
       $(SP_LIBd) \
       $(SIGIO_LIB4) \
       $(SFCIO_LIB4)

# If we can build wgrib2 on theia so that its static library does not 
# include NetCDF, then we can include WGRIB2_LIB and WGRIB2API_LIB as 
# part of the LIBS variable.  Otherwise, they will have to be added se-
# parately at the end of the build command below (as is currently done).
#        $(WGRIB2API_LIB) \
#        $(WGRIB2_LIB)

#
# Locations to search for module files referenced in USE statements and
# include files referenced in INCLUDE statements.
#
INCS = -I$(NEMSIO_INC) \
       -I$(WGRIB2API_INC) \
       -I$(SFCIO_INC4) \
       -I$(SIGIO_INC4) 

CMD = chgres_cube.exe

OBJS = chgres.o \
       atmosphere.o \
       input_data.o \
       grib2_util.o \
       surface.o \
       model_grid.o \
       program_setup.o \
       write_data.o \
       search_util.o \
       static_data.o \
       utils.o \
       atmdata_type.o
#
# Fortran compiler and compilation flags.
#
FCOMP = mpiifort
FFLAGS = -O0 -g -traceback -r8 -i4 -qopenmp -convert big_endian \
         -check bounds -warn unused -assume byterecl
#ESMF_F90ESMFLINKLIBS = -lesmf -cxxlib -lrt -ldl -lnetcdff -lnetcdf
#ESMF_F90ESMFLINKLIBS = -lesmf -cxxlib -lrt -ldl -lnetcdff
# -lnetcdf


$(CMD):	$(OBJS)
	$(FCOMP) $(FFLAGS) \
	  $(ESMF_F90COMPILEPATHS) \
	  -o $(CMD) $(OBJS) \
	  $(LIBS) \
	  $(ESMF_F90LINKPATHS) \
	  $(ESMF_F90ESMFLINKRPATHS) \
	  $(ESMF_F90ESMFLINKLIBS) \
	  $(WGRIB2_LIB) \
	  $(WGRIB2API_LIB) 

# We may be able to remove the above two lines involving WGRIB2_LIB and
# WGRIB2API_LIB if these can be made part of the definition of the vari-
# able LIBS above.  That depends on if we can get wgrib2 without NetCDF
# dependency to be installed on theia.
#
#	  $(WGRIB2_LIB)                # This is needed by WGRIB2API_LIB


model_grid.o: program_setup.o model_grid.F90
	$(FCOMP) $(FFLAGS) $(INCS) $(ESMF_F90COMPILEPATHS) -c model_grid.F90

utils.o: utils.f90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c utils.f90

atmdata_type.o: atmdata_type.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c atmdata_type.F90

grib2_util.o: model_grid.o program_setup.o atmdata_type.o grib2_util.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c grib2_util.F90

program_setup.o: program_setup.f90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c program_setup.f90

atmosphere.o: program_setup.o model_grid.o input_data.o atmosphere.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c atmosphere.F90

chgres.o: atmosphere.o model_grid.o program_setup.o surface.o chgres.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c chgres.F90

write_data.o: atmosphere.o model_grid.o program_setup.o surface.o static_data.o write_data.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c write_data.F90

input_data.o: program_setup.o model_grid.o atmdata_type.o grib2_util.o input_data.F90
	$(FCOMP) $(FFLAGS) $(INCS) $(ESMF_F90COMPILEPATHS) -c input_data.F90

surface.o: search_util.o model_grid.o input_data.o program_setup.o static_data.o surface.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c surface.F90

search_util.o: search_util.f90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c search_util.f90

static_data.o: model_grid.o program_setup.o static_data.F90
	$(FCOMP) $(FFLAGS) $(ESMF_F90COMPILEPATHS) -c static_data.F90

clean:
	rm -f *.o *.mod $(CMD) ../../exec/$(CMD)

install:
	mv $(CMD) ../../exec
