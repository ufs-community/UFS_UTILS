#!/bin/ksh

res=${res:-$1}
outdir=${outdir:-$2}

exec_dir=../exec
executable=$exec_dir/lakefrac
orodata_dir=/scratch1/NCEPDEV/nems/emc.nemspara/RT/NEMSfv3gfs/develop-20200603/INTEL

if [ $# -ne 2 ]; then 
echo "Usage $0 [res] [outdir]"
echo "res: 48, 96, ... , 1152, or all to create all resolutions."
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
  for res in 48 96 192 384 768 1152 
  do
    cp ${orodata_dir}/FV3_input_data_c${res}/INPUT/oro_data.tile?.nc .
    lakefrac 0 ${res}
  done
elif [ $res -eq 96 ]; then
  cp ${orodata_dir}/FV3_input_data/INPUT/oro_data.tile?.nc .
  lakefrac 0 ${res}
else
  cp ${orodata_dir}/FV3_input_data_c${res}/INPUT/oro_data.tile?.nc .
  lakefrac 0 ${res}
fi
