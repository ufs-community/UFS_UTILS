#!/bin/bash

#-----------------------------------------------------------------------------
# Invoke chgres to create 3-km CONUS coldstart files using HRRR GRIB2 data
# as input and with several surface variables from file instead of climo.
# The coldstart files are then compared to baseline files  using the 'nccmp' 
# utility.  This script is run by the machine specific driver script.
#-----------------------------------------------------------------------------

set -x

export DATA=$OUTDIR/3km_conus_hrrr_newsfc_grib2
rm -fr $DATA

export CRES=3357
export KMRES=3km
export FIXfv3=${HOMEreg}/fix/RRFS_CONUS_${KMRES}
export FIXsfc=${FIXfv3}/fix_sfc
export COMIN=${HOMEreg}/input_data/hrrr.grib2

export GRIB2_FILE_INPUT=1921300000000
export VCOORD_FILE=${HOMEufs}/fix/am/global_hyblev.l64.txt
export VARMAP_FILE=${HOMEufs}/parm/varmap_tables/GSDphys_var_map.txt
export INPUT_TYPE='grib2'
export CONVERT_NST=".false."
export OROG_FILES_TARGET_GRID="C3357_oro_data.tile7.nc"
export REGIONAL=1
export HALO_BLEND=0
export HALO_BNDY=4 
export CDATE=2019080100
export EXTERNAL_MODEL="HRRR"
export NSOILL_OUT=9
export TRACERS_TARGET='"NULL"'
export TRACERS_INPUT='"NULL"'
export SOTYP_FROM_CLIMO=.false.
export VGTYP_FROM_CLIMO=.false.
export VGFRC_FROM_CLIMO=.false.
export MINMAX_VGFRC_FROM_CLIMO=.false.
export TG3_FROM_SOIL=.true.
export LAI_FROM_CLIMO=.false.
export GEOGRID_FILE_INPUT=${HOMEufs}/fix/am/geo_em.d01.nc_HRRRX

export OMP_NUM_THREADS_CH=${OMP_NUM_THREADS:-1}

#-----------------------------------------------------------------------------
# Invoke chgres program.
#-----------------------------------------------------------------------------

echo "Starting at: " `date`

${HOMEufs}/ush/chgres_cube.sh

iret=$?
if [ $iret -ne 0 ]; then
  set +x
  echo "<<< 3-km CONUS HRRR W/ GSD PHYSICS AND SFC FROM FILE GRIB2 TEST FAILED. <<<"
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
    $NCCMP -dmfqS $files $HOMEreg/baseline_data/3km_conus_hrrr_newsfc_grib2/$files
    iret=$?
    if [ $iret -ne 0 ]; then
      test_failed=1
    fi
  fi
done

set +x
if [ $test_failed -ne 0 ]; then
  echo "<<< 3-km CONUS HRRR W/ GSD PHYSICS AND SFC FROM FILE GRIB2 TEST FAILED. >>>"
  if [ "$UPDATE_BASELINE" = "TRUE" ]; then
    $HOMEufs/reg_tests/update_baseline.sh $HOMEreg "3km_conus_hrrr_newsfc_grib2" $commit_num
  fi
else
  echo "<<< 3-km CONUS HRRR W/ GSD PHYSICS AND SFC FROM FILE GRIB2 TEST PASSED. >>>"
fi

exit 0
