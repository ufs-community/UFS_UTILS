#!/bin/ksh

#-----------------------------------------------------------
# MUST BE RUN WITH A MULTIPLE OF SIX MPI TASKS.  
# For regional grids, use any number of tasks.
#-----------------------------------------------------------

set -x

ulimit -s unlimited
ulimit -a

machine=${machine:-WCOSS_C}
res=${res:-96}
WORK_DIR=${WORK_DIR:-/scratch3/NCEPDEV/stmp1/$LOGNAME/gridgen_sfc.C${res}}
SAVE_DIR=${SAVE_DIR:-$WORK_DIR}
BASE_DIR=${BASE_DIR:-/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs}
EXEC_DIR=${EXEC_DIR:-$BASE_DIR/exec}
GRIDTYPE=${GRIDTYPE:-NULL}
FIX_FV3=${FIX_FV3:-/scratch4/NCEPDEV/global/save/glopara/svn/fv3gfs/fix_fv3/C$res}
mosaic_file=${mosaic_file:-$FIX_FV3/C${res}_mosaic.nc}
HALO=${HALO:-0}

if [ ! -d $SAVE_DIR ]; then
  mkdir -p $SAVE_DIR
fi

rm -fr $WORK_DIR
mkdir -p $WORK_DIR
cd $WORK_DIR

if [[ $GRIDTYPE == "nest" ]] || [[ $GRIDTYPE == "regional" ]]; then
  the_orog_files='"C'${res}'_oro_data.tile7.nc"'
else
  the_orog_files='"C'${res}'_oro_data.tile1.nc","C'${res}'_oro_data.tile2.nc","C'${res}'_oro_data.tile3.nc","C'${res}'_oro_data.tile4.nc","C'${res}'_oro_data.tile5.nc","C'${res}'_oro_data.tile6.nc"'
fi

if [[ $machine == "theia" ]]; then
  input_dir=/scratch4/NCEPDEV/da/noscrub/George.Gayno/climo_fields_netcdf
else
  input_dir=/gpfs/dell2/emc/modeling/noscrub/George.Gayno/landutil.git/climo_fields_netcdf
fi

cat << EOF > ./fort.41
&config
input_facsf_file="${input_dir}/facsf.1.0.nc"
input_substrate_temperature_file="${input_dir}/substrate_temperature.1.0.nc"
input_maximum_snow_albedo_file="${input_dir}/maximum_snow_albedo.0.05.nc"
input_snowfree_albedo_file="${input_dir}/snowfree_albedo.4comp.0.05.nc"
input_slope_type_file="${input_dir}/slope_type.1.0.nc"
input_soil_type_file="${input_dir}/soil_type.statsgo.0.05.nc"
input_vegetation_type_file="${input_dir}/vegetation_type.igbp.0.05.nc"
input_vegetation_greenness_file="${input_dir}/vegetation_greenness.0.144.nc"
mosaic_file_mdl="$mosaic_file"
orog_dir_mdl="$FIX_FV3"
orog_files_mdl=$the_orog_files
halo=$HALO
maximum_snow_albedo_method="bilinear"
snowfree_albedo_method="bilinear"
vegetation_greenness_method="bilinear"
/
EOF

if [[ $machine == "theia" ]]; then
  np=$PBS_NP
  module purge
  module load intel/15.1.133
  module load impi/5.1.1.109 
  module load netcdf/4.3.0
  APRUN_SFC=$(APRUN_SFC:-"mpirun -np $np")
  $APRUN_SFC $EXEC_DIR/gridgen_sfc.fv3
elif [[ $machine == "WCOSS_C" ]]; then
  APRUN_SFC=${APRUN_SFC:-"aprun -j 1 -n 6 -N 6"}
  $APRUN_SFC $EXEC_DIR/gridgen_sfc.fv3
fi

rc=$?

if [[ $rc == 0 ]]; then
  if [[ $GRIDTYPE != "regional" ]]; then
    for files in *.nc
    do
      if [[ -f $files ]]; then
        mv $files ${SAVE_DIR}/C${res}.${files}
      fi
    done
  else
    for files in *.halo.nc
    do
      if [[ -f $files ]]; then
        file2=${files%.halo.nc}
        mv $files ${SAVE_DIR}/C${res}.${file2}.halo${HALO}.nc
      fi
    done
    for files in *.nc
    do
      if [[ -f $files ]]; then
        file2=${files%.nc}
        mv $files ${SAVE_DIR}/C${res}.${file2}.halo0.nc
      fi
    done
  fi  # is regional?
else
  exit $rc
fi

exit 0
