#!/bin/bash
set -ex

#--Make symbolic links to 'fixed' directories.

RUN_ENVIR=${1}
machine=${2}

if [ $# -lt 2 ]; then
    set +x
    echo '***ERROR*** must specify two arguements: (1) RUN_ENVIR, (2) machine'
    echo ' Syntax: link_fv3gfs.sh ( nco | emc ) ( wcoss2 |  hera  | jet | orion | s4 )'
    exit 1
fi

if [ $RUN_ENVIR != emc -a $RUN_ENVIR != nco ]; then
    set +x
    echo '***ERROR*** unsupported run environment'
fi

if [ $machine != wcoss2 -a $machine != hera -a $machine != jet -a $machine != orion -a $machine != s4 ]; then
    set +x
    echo '***ERROR*** unsupported machine'
    echo 'Syntax: link_fv3gfs.sh ( nco | emc ) ( wcoss2 | hera | jet | orion | s4 )'
    exit 1
fi

LINK="ln -fs"
SLINK="ln -fs"
[[ $RUN_ENVIR = nco ]] && LINK="cp -rp"

pwd=$(pwd -P)

#------------------------------
#--model fix fields
#------------------------------
if [ $machine = "hera" ]; then
    FIX_DIR="/scratch1/NCEPDEV/global/glopara/fix"
elif [ $machine = "jet" ]; then
    FIX_DIR="/lfs4/HFIP/hfv3gfs/glopara/git/fv3gfs/fix"
elif [ $machine = "orion" ]; then
    FIX_DIR="/work/noaa/global/glopara/fix"
elif [ $machine = "wcoss2" ]; then
    FIX_DIR="/lfs/h2/emc/global/save/emc.global/FIX/fix"
elif [ $machine = "s4" ]; then
    FIX_DIR="/data/prod/glopara/fix"
fi

am_ver=${am_ver:-20220805}
orog_ver=${orog_ver:-20220805}
sfc_climo_ver=${sfc_climo_ver:-20220805}

for dir in am orog sfc_climo; do
    if [ -d $dir ]; then
      [[ $RUN_ENVIR = nco ]] && chmod -R 755 $dir
      rm -rf $dir
    fi
    fix_ver="${dir}_ver"
    $LINK $FIX_DIR/$dir/${!fix_ver} ${dir}
done

exit 0
