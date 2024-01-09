#!/bin/bash

#-----------------------------------------------------------------------------
# Invoke chgres to create 25-km CONUS coldstart files using GFS 
# PGRIB2+BGRIB2 data as input.  The coldstart files are then compared to
# baseline files using the 'nccmp' utility.  This script is run by the
# This script is run by the machine specific driver script.
#-----------------------------------------------------------------------------

set -x

export DATA=$OUTDIR/25km_conus_gfs_pbgrib2
rm -fr $DATA

export CRES=405
export KMRES=25km
export FIXfv3=${HOMEreg}/fix.new/RRFS_CONUS_${KMRES}

export COMIN=${HOMEreg}/input_data/gfs.pbgrib2

export GRIB2_FILE_INPUT=gfs.t18z.pgrb2.0p25.f006
export VCOORD_FILE=${HOMEufs}/fix/am/global_hyblev.l64.txt
export VARMAP_FILE=${HOMEufs}/parm/varmap_tables/GFSphys_var_map.txt
export INPUT_TYPE='grib2'
export CONVERT_NST=".false."
export OROG_FILES_TARGET_GRID="C405_oro_data.tile7.nc"
export REGIONAL=1
export HALO_BLEND=0
export HALO_BNDY=4 
export CDATE=2021062718

export OMP_NUM_THREADS_CH=${OMP_NUM_THREADS:-1}

#-----------------------------------------------------------------------------
# Invoke chgres program.
#-----------------------------------------------------------------------------

echo "Starting at: " `date`

${HOMEufs}/ush/chgres_cube.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< 25-km CONUS GFS PGRIB2+BGRIB2 TEST FAILED. >>>"
  exit $iret
fi

echo "Ending at: " `date`

#-----------------------------------------------------------------------------
# Compare output from chgres to baseline set of data.
#
# orion's nccmp utility does not work with the netcdf
# required to run ufs_utils.  So swap it.
#-----------------------------------------------------------------------------

machine=${machine:-NULL}
if [ $machine == 'orion' ]; then
  module unload netcdfp/4.7.4.release
  module load netcdf/4.7.2
fi

cd $DATA

test_failed=0
for files in *.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/25km_conus_gfs_pbgrib2/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< 25-KM CONUS GFS PGRIB2+BGRIB2 TEST FAILED. >>>"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    $HOMEufs/reg_tests/update_baseline.sh $HOMEreg "25km_conus_gfs_pbgrib2" $commit_num
  fi
else
  echo "<<< 25-KM CONUS GFS PGRIB2+BGRIB2 TEST PASSED. >>>"
fi

exit 0
