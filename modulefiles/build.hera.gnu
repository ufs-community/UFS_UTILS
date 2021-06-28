#%Module#####################################################
## Build and run module for Hera
#############################################################

module load hpss
module load cmake/3.16.1

module use /scratch2/NCEPDEV/nwprod/hpc-stack/libs/hpc-stack/modulefiles/stack

module load hpc/1.1.0
module load hpc-gnu/9.2.0
module load hpc-mpich/3.3.2

module load netcdf/4.7.4
module load esmf/8_1_0_beta_snapshot_27
module load bacio/2.4.1
module load g2/3.4.1
module load ip/3.3.3
module load nemsio/2.5.2
module load sp/2.3.3
module load w3nco/2.4.1
module load sfcio/1.4.1
module load sigio/2.3.2
module load wgrib2/2.0.8
module load nccmp/1.8.7.0
module load png/1.6.35
module load zlib/1.2.11
module load jasper/2.0.22
