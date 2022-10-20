EXEC    =soil.exe
FC	=ftn
FOPTS	= -O0 -FR -I${NETCDF}/include -assume byterecl -warn unused -check bounds
OBJ_TST =soil.o
OBJS=${OBJ_TST}
# *************************************************************************
all: chop
                                                                                
chop : ${OBJS}
	$(FC) ${FOPTS} -o ${EXEC} -L${NETCDF}/lib -lnetcdf -lnetcdff -L${HDF5_ROOT}/lib -lhdf5 -lhdf5_fortran ${OBJS}

                                                                                
.SUFFIXES : .o .f90

.f90.o :
	$(FC) ${FOPTS} -c $*.f90

clean:
	rm -f *.o *.mod $(EXEC)
