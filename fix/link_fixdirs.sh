#!/bin/bash
set -ex

# Set up the 'fixed' directories. 
# 
# This script takes two arguments:
#
#  $RUN_ENVIR - Either 'emc' (creates links) or
#               'nco' (copies data).
#
#  $machine - is the machine. Choices are:
#             'wcoss2', 'hera', 'jet', 'orion', 'hercules', 's4', 'gaea'

RUN_ENVIR=${1}
machine=${2}

if [ $# -lt 2 ]; then
    set +x
    echo '***ERROR*** must specify two arguements: (1) RUN_ENVIR, (2) machine'
    echo ' Syntax: link_fv3gfs.sh ( nco | emc ) ( wcoss2 |  hera  | jet | orion | hercules | s4 | gaea )'
    exit 1
fi

if [ $RUN_ENVIR != emc -a $RUN_ENVIR != nco ]; then
    set +x
    echo '***ERROR*** unsupported run environment'
    echo ' Must choose either "nco" or "emc".'
    exit 1
fi

if [ $machine != wcoss2 -a $machine != hera -a $machine != jet -a $machine != orion -a $machine != s4 -a $machine != hercules -a $machine != gaea ]; then
    set +x
    echo '***ERROR*** unsupported machine'
    echo 'Syntax: link_fv3gfs.sh ( nco | emc ) ( wcoss2 | hera | jet | orion | hercules | s4 | gaea )'
    exit 1
fi

LINK="ln -fs"
SLINK="ln -fs"
[[ $RUN_ENVIR = nco ]] && LINK="cp -rpL"

pwd=$(pwd -P)

#------------------------------
#--model fix fields
#------------------------------
if [ $machine = "hera" ]; then
    FIX_DIR="/scratch1/NCEPDEV/global/glopara/fix"
elif [ $machine = "jet" ]; then
    FIX_DIR="/lfs4/HFIP/hfv3gfs/glopara/git/fv3gfs/fix"
elif [ $machine = "orion" -o $machine = "hercules" ]; then
    FIX_DIR="/work/noaa/global/glopara/fix"
elif [ $machine = "wcoss2" ]; then
    FIX_DIR="/lfs/h2/emc/global/noscrub/emc.global/FIX/fix"
elif [ $machine = "s4" ]; then
    FIX_DIR="/data/prod/glopara/fix"
elif [ $machine = "gaea" ]; then
    FIX_DIR="/gpfs/f5/epic/proj-shared/global/glopara/data/fix"
fi

am_ver=${am_ver:-20220805}
orog_ver=${orog_ver:-20231027}
sfc_climo_ver=${sfc_climo_ver:-20230925}

for dir in am orog orog_raw sfc_climo; do
    if [ -d $dir ]; then
      [[ $RUN_ENVIR = nco ]] && chmod -R 755 $dir
      rm -rf $dir
    fi
    if [ $dir = "orog_raw" ]; then
      $LINK $FIX_DIR/raw/orog ${dir}
    else
      fix_ver="${dir}_ver"
      $LINK $FIX_DIR/$dir/${!fix_ver} ${dir}
    fi
done

exit 0
