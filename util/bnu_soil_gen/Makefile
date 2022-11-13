#

NETCDF_INC = -I${NETCDF}/include
NETCDF_LIB = -L${NETCDF}/lib -lnetcdf -lnetcdff

F90 = ifort

soil_fill: soil_fill.f90
	$(F90) $(FFLAGS) -o $(@) $(NETCDF_INC) soil_fill.f90 $(NETCDF_LIB)

