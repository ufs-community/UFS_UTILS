#!/bin/bash

#-----------------------------------------------------------------------------
# Invoke chgres to create C192 coldstart files using GFS GRIB2 data
# as input.  The coldstart files are then compared to baseline files
# using the 'nccmp' utility.  This script is run by the machine specific 
# driver script.
#-----------------------------------------------------------------------------

set -x

export DATA=$OUTDIR/c192_gfs_grib2
rm -fr $DATA

export CRES=192
export FIXfv3=${HOMEreg}/fix/C192
export COMIN=${HOMEreg}/input_data/gfs.grib2

export GRIB2_FILE_INPUT=gfs.t00z.pgrb2.0p50.f000
export VCOORD_FILE=${HOMEufs}/fix/fix_am/global_hyblev.l65.txt
export VARMAP_FILE=${HOMEufs}/parm/varmap_tables/GFSphys_var_map.txt
export INPUT_TYPE='grib2'
export CONVERT_NST=".false."

export CDATE=2019110400

export OMP_NUM_THREADS_CH=${OMP_NUM_THREADS:-1}

NCCMP=${NCCMP:-$(which nccmp)}

#-----------------------------------------------------------------------------
# Invoke chgres program.
#-----------------------------------------------------------------------------

echo "Starting at: " `date`

${HOMEufs}/ush/chgres_cube.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< C192 GFS GRIB2 TEST FAILED. <<<"
  exit $iret
fi

echo "Ending at: " `date`

#-----------------------------------------------------------------------------
# Compare output from chgres to baseline set of data.
#-----------------------------------------------------------------------------

cd $DATA

test_failed=0
for files in *.nc
do
  if [ -f $files ]; then
    echo CHECK $files
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/c192_gfs_grib2/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< C192 GFS GRIB2 TEST FAILED. >>>"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    $HOMEufs/reg_tests/update_baseline.sh $HOMEreg "c192_gfs_grib2" $commit_num
  fi
else
  echo "<<< C192 GFS GRIB2 TEST PASSED. >>>"
fi

exit 0
