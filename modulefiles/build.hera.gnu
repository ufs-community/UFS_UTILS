#%Module#####################################################
## Build and run module for Hera
#############################################################

module load hpss
module load cmake/3.16.1
module load gnu/9.2.0
module use -a /scratch1/BMC/gmtb/software/modulefiles/gnu-9.2.0/mpich-3.3.2
module load mpich/3.3.2

module use -a /scratch1/BMC/gmtb/software/ufs-stack-20200909/gnu-9.2.0/mpich-3.3.2/modules

module load libpng/1.6.35
module load netcdf/4.7.4
module load esmf/8.1.0bs27

module load bacio/2.4.0
module load g2/3.4.0
module load ip/3.3.0
module load nemsio/2.5.1
module load sp/2.3.0
module load w3emc/2.7.0
module load w3nco/2.4.0
module load gfsio/1.4.0
module load sfcio/1.4.0
module load sigio/2.3.0
module load nemsiogfs/2.5.0
module load landsfcutil/2.4.0
module load wgrib2/2.0.8
