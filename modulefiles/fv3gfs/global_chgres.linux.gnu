#############################################################
## global_chgres component - linux.gnu
#############################################################

export GFSIO_INC4=${NCEPLIBS}/gfsio/include_4
export IP_INCd=${NCEPLIBS}/ip/include_d
export LANDSFCUTIL_INCd=${NCEPLIBS}/landsfcutil/include
export NEMSIOGFS_INC=${NCEPLIBS}/nemsiogfs/include
export NEMSIO_INC=${NCEPLIBS}/nemsio/include
export SFCIO_INC4=${NCEPLIBS}/sfcio/include
export SIGIO_INC4=${NCEPLIBS}/sigio/include

export BACIO_LIB4=${NCEPLIBS}/bacio/lib/libbacio_v2.1.0_4.a
export GFSIO_LIB4=${NCEPLIBS}/gfsio/lib/libgfsio_v1.1.0_4.a
export IP_LIBd=${NCEPLIBS}/ip/lib/libip_v3.0.0_d.a
export LANDSFCUTIL_LIBd=${NCEPLIBS}/landsfcutil/lib/liblandsfcutil_v2.1.0_d.a
export NEMSIOGFS_LIB=${NCEPLIBS}/nemsiogfs/lib/libnemsiogfs_v2.2.1.a
export NEMSIO_LIB=${NCEPLIBS}/nemsio/lib/libnemsio_v2.2.3.a
export SFCIO_LIB4=${NCEPLIBS}/sfcio/lib/libsfcio_v1.1.0_4.a
export SIGIO_LIB4=${NCEPLIBS}/sigio/lib/libsigio_v2.1.0_4.a
export SP_LIBd=${NCEPLIBS}/sp/lib/libsp_v2.0.2_d.a
export W3EMC_LIBd=${NCEPLIBS}/w3emc/lib/libw3emc_v2.3.0_d.a
export W3NCO_LIBd=${NCEPLIBS}/w3nco/lib/libw3nco_v2.0.6_d.a

export NETCDF_INCLUDE="-I${NETCDF}/include"
export NETCDF_LDFLAGS_F="-L${NETCDF}/lib -lnetcdff -lnetcdf -L${HDF5}/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl -lz -lm"

export FCMP=gfortran
export FFLAGSM="-O3 -fbacktrace -fdefault-real-8 -fconvert=big-endian"
export LDFLAGSM="-fopenmp"
export OMPFLAGM="-fopenmp"
