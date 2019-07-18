#%Module#####################################################
## chgres build module for Odin
#############################################################

module use /oldscratch/ywang/external/modulefiles
module load esmf/8.0.0bs30

module load cray-netcdf-hdf5parallel
module load cray-parallel-netcdf
module load cray-hdf5-parallel
module load w3nco/v2.0.6
module load nemsio/v2.2.2
module load bacio/v2.0.2
module load sp/v2.0.2
module load sfcio/v1.0.0
module load sigio/v2.0.1


export FCOMP=ftn
export FFLAGS="-O3 -fp-model precise -g -traceback -r8 -i4 -qopenmp -convert big_endian -assume byterecl"
export WGRIB2API_LIB="/home/larissa.reames/tmp/wgrib2-2/grib2/lib/libwgrib2_api.a"
export WGRIB2API_INC="/home/larissa.reames/tmp/wgrib2-2/grib2/lib"
export WGRIB2_LIB="/home/larissa.reames/tmp/wgrib2-2/grib2/lib/libwgrib2.a"


# for debugging
#export FFLAGS="-O0 -g -traceback -r8 -i4 -qopenmp -convert big_endian -check bounds -warn unused -assume byterecl"
