#############################################################
## chgres_cube component - linux.gnu
#############################################################

export IP_INCd=${NCEPLIBS}/ip/include_d
export NEMSIO_INC=${NCEPLIBS}/nemsio/include
export SFCIO_INC4=${NCEPLIBS}/sfcio/include
export SIGIO_INC4=${NCEPLIBS}/sigio/include

export BACIO_LIB4=${NCEPLIBS}/bacio/lib/libbacio_v2.1.0_4.a
export IP_LIBd=${NCEPLIBS}/ip/lib/libip_v3.0.0_d.a
export NEMSIO_LIB=${NCEPLIBS}/nemsio/lib/libnemsio_v2.2.3.a
export SFCIO_LIB4=${NCEPLIBS}/sfcio/lib/libsfcio_v1.1.0.a
export SIGIO_LIB4=${NCEPLIBS}/sigio/lib/libsigio_v2.1.0.a
export SP_LIBd=${NCEPLIBS}/sp/lib/libsp_v2.0.2_d.a
export W3NCO_LIBd=${NCEPLIBS}/w3nco/lib/libw3nco_v2.0.6_d.a

export FCOMP=mpif90
export FFLAGS="-O3 -g -fbacktrace -fdefault-real-8 -ffree-line-length-none -fopenmp -fconvert=big-endian"
# for debugging
#export FFLAGS="-O0 -g -fbacktrace -fdefault-real-8 -ffree-line-length-none -fopenmp -fconvert=big-endian"
