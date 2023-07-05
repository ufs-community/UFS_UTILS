#!/bin/bash

#
# Shan's original code  modified to accept arguments 
# Land points will change hence orog is re run





nargv=$#

if [ $nargv -eq 2 ];  then  # we have the right number of variables

 CRES=$1 
 ocn=$2

else
  echo "Number of arguments must be 2 for fv3gfs_ocean_merge"
  echo "orog.sh res ocn"
  exit

fi

UFSDIR=$home_dir


source $UFSDIR/sorc/machine-setup.sh > /dev/null 2>&1
module use $UFSDIR/modulefiles
module load build.hera.intel
module list

FIXDIR=$UFSDIR/fix

executable=$UFSDIR/exec/orog
exec_filter=$UFSDIR/exec/filter_topo

  
SHAN_DIR=$TEMP_DIR/ocean_merged/C${CRES}.mx${ocn}
WORKDIR=$(dirname $out_dir)/C${CRES}.mx${ocn}/

rm -rf $WORKDIR
mkdir -p $WORKDIR
cd $WORKDIR



ln -fs $FIXDIR/orog/thirty.second.antarctic.new.bin fort.15
ln -fs $FIXDIR/orog/landcover30.fixed .
ln -fs $FIXDIR/orog/gmted2010.30sec.int fort.235

for tiles in tile1 tile2 tile3 tile4 tile5 tile6
do
#ln -fs $FIXDIR/fix_fv3_gmted2010/C${CRES}/C${CRES}_grid.${tiles}.nc .

ln -fs $FIXDIR/orog/C${CRES}/C${CRES}_grid.${tiles}.nc .
cp ${SHAN_DIR}/C${CRES}.mx${ocn}.${tiles}.nc ${SHAN_DIR}/oro_data.${tiles}.nc



cat << EOF >./fort.41
 &test
   mom6_file="${SHAN_DIR}/oro_data.${tiles}.nc"
 /  
EOF

cat << EOF >./INPS
1 ${CRES} ${CRES} 0 0 0 0 0 0
C${CRES}_grid.${tiles}.nc
none
EOF

cat INPS
time $executable < INPS

iret=$?
if [ $iret != 0 ]; then
  echo error
  exit $iret
fi

#cp out.oro.nc  oro.C${CRES}.${tiles}.no.filter.nc
mv out.oro.nc  C${CRES}_oro_data.${tiles}.nc

rm -f ./fort.41
rm -f ./INPS

done

# run filtering

cat << EOF >./input.nml
 &filter_topo_nml
   grid_file="$FIXDIR/orog/C${CRES}/C${CRES}_mosaic.nc"
   topo_file="C${CRES}_oro_data"
   mask_field="land_frac"
   regional=.false.
   stretch_fac=1.0
   res=$CRES
 /  
EOF

time ${exec_filter}

iret=$?
if [ $iret != 0 ]; then
  echo error
  exit $iret
fi

cp $out_dir/C${CRES}_mosaic.nc $WORKDIR 

