#############################################################
## global_chgres component - linux.gnu
#############################################################

export SIGIO_INC4=${NCEPLIBS}/sigio/include/sigio_v2.1.0_4
export SFCIO_INC4=${NCEPLIBS}/sfcio/include/sfcio_v1.1.0_4
export LANDSFCUTIL_INCd=${NCEPLIBS}/landsfcutil/include/landsfcutil_v2.1.0_d
export NEMSIO_INC=${NCEPLIBS}/nemsio/include/nemsio_v2.2.5
export NEMSIOGFS_INC=${NCEPLIBS}/nemsiogfs/include/nemsiogfs_v2.2.0
export GFSIO_INC4=${NCEPLIBS}/gfsio/include/gfsio_v1.1.0_4
export IP_INCd=${NCEPLIBS}/ip/include/ip_v3.0.1_d

export GFSIO_LIB4=${NCEPLIBS}/gfsio/libgfsio_v1.1.0_4.a
export NEMSIOGFS_LIB=${NCEPLIBS}/nemsiogfs/libnemsiogfs_v2.2.0.a
export NEMSIO_LIB=${NCEPLIBS}/nemsio/libnemsio_v2.2.5.a
export SIGIO_LIB4=${NCEPLIBS}/sigio/libsigio_v2.1.0_4.a
export SFCIO_LIB4=${NCEPLIBS}/sfcio/libsfcio_v1.1.0_4.a
export LANDSFCUTIL_LIBd=${NCEPLIBS}/landsfcutil/liblandsfcutil_v2.1.0_d.a
export IP_LIBd=${NCEPLIBS}/ip/libip_v3.0.1_d.a
export SP_LIBd=${NCEPLIBS}/sp/libsp_v2.0.2_d.a
export W3EMC_LIBd=${NCEPLIBS}/w3emc/libw3emc_v2.3.0_d.a
export W3NCO_LIBd=${NCEPLIBS}/w3nco/libw3nco_v2.0.6_d.a
export BACIO_LIB4=${NCEPLIBS}/bacio/libbacio_v2.0.2_4.a

export NETCDF_INCLUDE="-I${NETCDF}/include"
export NETCDF_LDFLAGS_F="-L${NETCDF}/lib -lnetcdff -lnetcdf -L${HDF5}/lib -lhdf5_fortran -lhdf5_hl -lhdf5 -ldl -lz -lm"

export FCMP=gfortran
export FFLAGSM="-O3 -fbacktrace -fdefault-real-8 -fconvert=big-endian"
export LDFLAGSM="-fopenmp"
export OMPFLAGM="-fopenmp"
