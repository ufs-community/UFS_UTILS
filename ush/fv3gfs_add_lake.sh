#!/bin/ksh

res=${res:-$1}
outdir=${outdir:-$2}

exec_dir=../exec
executable=$exec_dir/lakefrac
orodata_dir=/scratch1/NCEPDEV/nems/emc.nemspara/RT/NEMSfv3gfs/develop-20200603/INTEL

if [ $# -ne 2 ]; then 
echo "Usage $0 [res] [outdir]"
echo "res: 96, ... , 768, or all to create all resolutions."
echo "outdir: output directory."
exit 1 
fi

if [ ! -s $executable ]; then
  set +x
  echo
  echo "FATAL ERROR: ${executable} does not exist."
  echo
  set -x
  exit 1 
fi

if [ ! -s $outdir ]; then  mkdir -p $outdir ;fi

cp ${executable}  ${outdir}
cd $outdir

if [ $res -eq "all" ]; then 
  for res in 96 192 384 768
  do
    if [ $res -eq 96 ]; then 
      cp ${orodata_dir}/FV3_input_data/INPUT/oro_data.tile?.nc .
    else
      cp ${orodata_dir}/FV3_input_data_c${res}/INPUT/oro_data.tile?.nc .
    fi
    lakefrac 0 ${res}
    for tile in 1 2 3 4 5 6
    do
      mv oro_data.tile${tile}.nc C${res}_oro_data.tile${tile}.nc
    done 
  done
elif [ $res -eq 96 ]; then
  cp ${orodata_dir}/FV3_input_data/INPUT/oro_data.tile?.nc .
  lakefrac 0 ${res}
  for tile in 1 2 3 4 5 6
  do
    mv oro_data.tile${tile}.nc C${res}_oro_data.tile${tile}.nc
  done 
else
  cp ${orodata_dir}/FV3_input_data_c${res}/INPUT/oro_data.tile?.nc .
  lakefrac 0 ${res}
  for tile in 1 2 3 4 5 6
  do
    mv oro_data.tile${tile}.nc C${res}_oro_data.tile${tile}.nc
  done 
fi
