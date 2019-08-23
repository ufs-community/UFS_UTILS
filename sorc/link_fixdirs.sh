#!/bin/ksh
set -ex

#--Make symbolic links to 'fixed' directories.

RUN_ENVIR=${1}
machine=${2}

if [ $# -lt 2 ]; then
    set +x
    echo '***ERROR*** must specify two arguements: (1) RUN_ENVIR, (2) machine'
    echo ' Syntax: link_fv3gfs.sh ( nco | emc ) ( cray | dell | theia | jet )'
    exit 1
fi

if [ $RUN_ENVIR != emc -a $RUN_ENVIR != nco ]; then
    set +x
    echo '***ERROR*** unsupported run environment'
    echo 'Syntax: link_fv3gfs.sh ( nco | emc ) ( cray | dell | theia | jet )'
    exit 1
fi
if [ $machine != cray -a $machine != theia -a $machine != dell -a $machine != jet ]; then
    set +x
    echo '***ERROR*** unsupported machine'
    echo 'Syntax: link_fv3gfs.sh ( nco | emc ) ( cray | dell | theia | jet )'
    exit 1
fi

LINK="ln -fs"
SLINK="ln -fs"
[[ $RUN_ENVIR = nco ]] && LINK="cp -rp"

pwd=$(pwd -P)

#------------------------------
#--model fix fields
#------------------------------
if [ $machine == "cray" ]; then
    FIX_DIR="/gpfs/hps3/emc/global/noscrub/emc.glopara/git/fv3gfs/fix"
elif [ $machine = "dell" ]; then
    FIX_DIR="/gpfs/dell2/emc/modeling/noscrub/emc.glopara/git/fv3gfs/fix"
elif [ $machine = "theia" ]; then
    FIX_DIR="/scratch4/NCEPDEV/global/save/glopara/git/fv3gfs/fix"
elif [ $machine = "jet" ]; then
    FIX_DIR="/lfs3/projects/hfv3gfs/glopara/git/fv3gfs/fix"
fi
cd ${pwd}/../fix                ||exit 8
for dir in fix_am fix_fv3 fix_orog fix_fv3_gmted2010 ; do
    [[ -d $dir ]] && rm -rf $dir
    $LINK $FIX_DIR/$dir  .
done

if [ $machine == "cray" ] || [ $machine = "dell" ]; then
    $LINK /gpfs/dell2/emc/modeling/noscrub/George.Gayno/landutil.git/climo_fields_netcdf ./fix_sfc_climo
elif [ $machine = "theia" ]; then
    $LINK /scratch4/NCEPDEV/da/noscrub/George.Gayno/climo_fields_netcdf ./fix_sfc_climo
elif [ $machine = "jet" ]; then
    $LINK /mnt/lfs3/projects/emcda/George.Gayno/climo_fields_netcdf ./fix_sfc_climo
fi

exit 0
