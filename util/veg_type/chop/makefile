EXEC    =vegt.exe
FC	=ifort
FOPTS	= -O0 -FR -I${NETCDF}/include -check bounds -assume noold_ldout_format -assume byterecl
OBJ_TST =vegt.o
OBJS=${OBJ_TST}
# *************************************************************************
all: vegt
                                                                                
vegt : ${OBJS}
	$(FC) ${FOPTS} -o ${EXEC} -L${NETCDF}/lib -lnetcdf -lnetcdff -L${HDF5}/lib -lhdf5 -lhdf5_fortran ${OBJS}

                                                                                
.SUFFIXES : .o .f90

.f90.o :
	$(FC) ${FOPTS} -c $*.f90

clean:
	rm -f *.o *.mod $(EXEC)
