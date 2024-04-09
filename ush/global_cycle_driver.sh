#!/bin/bash

set -eux

#-------------------------------------------------------------------------------------------------
# Update surface fields in FV3 restart files on the cubed-sphere grid
# George Gayno,  09/XX/2017
# Rahul Mahajan, 10/11/2017
#-------------------------------------------------------------------------------------------------

export CASE=${CASE:-C768}                    # resolution of tile: 48, 96, 192, 384, 768, 1152, 3072
export OCNRES=${OCNRES:-100}
export CDATE=${CDATE:-${cdate:-2017031900}}  # format yyyymmddhh yyyymmddhh ...
export CDUMP=${CDUMP:-gfs}                   # gfs or gdas
export COMPONENT=${COMPONENT:-atmos}

pwd=$(pwd)
export DMPDIR=${DMPDIR:-$pwd}
export PACKAGEROOT=${PACKAGEROOT:-/lfs/h1/ops/prod/packages}
export gfs_ver=${gfs_ver:-v15.0.0}
export HOMEgfs=${HOMEgfs:-${PACKAGEROOT}/gfs.${gfs_ver}}
export FIXgfs=${FIXgfs:-$HOMEgfs/fix}   

ntiles=${ntiles:-6}
DONST=${DONST:-"NO"}
COMIN=${COMIN:-$pwd}
COMOUT=${COMOUT:-$pwd}

CYCLESH=${CYCLESH:-$HOMEgfs/ush/global_cycle.sh}
export CYCLEXEC=${CYCLEXEC:-$HOMEgfs/exec/global_cycle}
export OMP_NUM_THREADS_CY=${OMP_NUM_THREADS_CY:-24}
export APRUNCY=${APRUNCY:-"time"}
export VERBOSE=${VERBOSE:-"YES"}

export FHOUR=${FHOUR:-0}
export DELTSFC=${DELTSFC:-6}

PDY=$(echo $CDATE | cut -c1-8)
cyc=$(echo $CDATE | cut -c9-10)

export FNTSFA=${FNTSFA:-$DMPDIR/${CDUMP}.${PDY}/${cyc}/${COMPONENT}/${CDUMP}.t${cyc}z.rtgssthr.grb}
export FNSNOA=${FNSNOA:-$DMPDIR/${CDUMP}.${PDY}/${cyc}/${COMPONENT}/${CDUMP}.t${cyc}z.snogrb_t1534.3072.1536}
export FNACNA=${FNACNA:-$DMPDIR/${CDUMP}.${PDY}/${cyc}/${COMPONENT}/${CDUMP}.t${cyc}z.seaice.5min.blend.grb}

export CYCLVARS=${CYCLVARS:-"FSNOL=-2.,FSNOS=99999.,"}

if [ $DONST = "YES" ]; then
    export NST_FILE=${NST_FILE:-$COMOUT/dtfanl.nc}
else
    export NST_FILE="NULL"
fi

export DO_SFCCYLE=${DO_SFCCYCLE:-".true."}
export DO_LNDINC=${DO_LNDINC:-".false."}
export DO_SOI_INC_GSI=${DO_SOI_INC_GSI:-".false."}
export DO_SNO_INC_JEDI=${DO_SNO_INC_JEDI:-".false."}
export DO_SOI_INC_JEDI=${DO_SOI_INC_JEDI:-".false."}
export FRAC_GRID=${FRAC_GRID:-".false."}

CRES=$(echo $CASE | cut -c 2-)
JCAP_CASE=$((2*CRES-2))
LONB_CASE=$((4*CRES))
LATB_CASE=$((2*CRES))

export JCAP=${JCAP:-$JCAP_CASE}
export LONB=${LONB:-$LONB_CASE}
export LATB=${LATB:-$LATB_CASE}

export MAX_TASKS_CY=${MAX_TASKS_CY:-99999}

# Temporary rundirectory
export DATA=${DATA:-$pwd/rundir$$}
rm -fr $DATA
mkdir -p $DATA

# Make a copy of the input restart file as fnbgso. For fractional grids,
# only a few fields will be updated by sfcsub.F. Only these records
# will be overwritten in fnbgso.

for n in $(seq 1 $ntiles); do
  ln -fs $COMIN/$PDY.${cyc}0000.sfc_data.tile${n}.nc      $DATA/fnbgsi.00$n

# Make a copy of the input restart file in the working directory.
# global_cycle will update the required records for noah-mp.

  cp $COMIN/$PDY.${cyc}0000.sfc_data.tile${n}.nc $COMOUT/$PDY.${cyc}0000.sfcanl_data.tile${n}.nc
  chmod 644  $COMOUT/$PDY.${cyc}0000.sfcanl_data.tile${n}.nc
  ln -fs $COMOUT/$PDY.${cyc}0000.sfcanl_data.tile${n}.nc  $DATA/fnbgso.00$n

  ln -fs $FIXgfs/orog/${CASE}/C${CRES}_grid.tile${n}.nc       $DATA/fngrid.00$n
  if (( OCNRES > 9999 ));then
    ln -fs $FIXgfs/orog/${CASE}/C${CRES}_oro_data.tile${n}.nc   $DATA/fnorog.00$n
  else
    ln -fs $FIXgfs/orog/${CASE}/C${CRES}.mx${OCNRES}_oro_data.tile${n}.nc   $DATA/fnorog.00$n
  fi

  if [[ "$DO_SNO_INC_JEDI" == ".true." ]] ; then  
        ln -fs $COMIN/$PDY.${cyc}0000.xainc.tile${n}.nc      $DATA/snow_xainc.00$n
  fi

  if [[ "$DO_SOI_INC_JEDI" == ".true." ]] ; then
        ln -fs $COMIN/soil_sfcincr_jedi.00$n $DATA/soil_xainc.00$n
  fi

  if [[ "$DO_SOI_INC_GSI" == ".true." ]] ; then
        ln -fs $COMIN/sfcincr_gsi.00$n $DATA/sfcincr_gsi.00$n
  fi
done

$CYCLESH

rc=$?
if [[ $rc -ne 0 ]] ; then
    echo "***ERROR*** rc= $rc"
    exit $rc
fi

exit 0
