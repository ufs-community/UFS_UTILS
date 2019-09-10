#############################################################
## chgres_cube component - linux.gnu
#############################################################

export SIGIO_INC4=${NCEPLIBS}/sigio/include/sigio_v2.1.0_4
export SFCIO_INC4=${NCEPLIBS}/sfcio/include/sfcio_v1.1.0_4
export NEMSIO_INC=${NCEPLIBS}/nemsio/include/nemsio_v2.2.5
export IP_INCd=${NCEPLIBS}/ip/include/ip_v3.0.1_d

export W3NCO_LIBd=${NCEPLIBS}/w3nco/libw3nco_v2.0.6_d.a
export BACIO_LIB4=${NCEPLIBS}/bacio/libbacio_v2.0.2_4.a
export NEMSIO_LIB=${NCEPLIBS}/nemsio/libnemsio_v2.2.5.a
export SIGIO_LIB4=${NCEPLIBS}/sigio/libsigio_v2.1.0_4.a
export SFCIO_LIB4=${NCEPLIBS}/sfcio/libsfcio_v1.1.0_4.a
export IP_LIBd=${NCEPLIBS}/ip/libip_v3.0.1_d.a
export SP_LIBd=${NCEPLIBS}/sp/libsp_v2.0.2_d.a

export FCOMP=mpif90
export FFLAGS="-O3 -g -fbacktrace -fdefault-real-8 -ffree-line-length-none -fopenmp -fconvert=big-endian"
# for debugging
#export FFLAGS="-O0 -g -fbacktrace -fdefault-real-8 -ffree-line-length-none -fopenmp -fconvert=big-endian"
