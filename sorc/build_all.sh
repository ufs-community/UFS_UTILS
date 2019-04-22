#!/bin/sh
set -eu
#------------------------------------
# USER DEFINED STUFF:
#
# USE_PREINST_LIBS: set to "true" to use preinstalled libraries.
#                   Anything other than "true"  will use libraries locally.
#------------------------------------

export USE_PREINST_LIBS="true"

#------------------------------------
# END USER DEFINED STUFF
#------------------------------------

build_dir=`pwd`
logs_dir=$build_dir/logs
if [ ! -d $logs_dir  ]; then
  echo "Creating logs folder"
  mkdir $logs_dir
fi

# Check final exec folder exists
if [ ! -d "../exec" ]; then
  echo "Creating ../exec folder"
  mkdir ../exec
fi

#------------------------------------
# INCLUDE PARTIAL BUILD 
#------------------------------------

. ./partial_build.sh

#------------------------------------
# build NEMS util
#------------------------------------
$Build_nems_util && {
echo " .... Building NEMS util .... "
./build_nems_util.sh > $logs_dir/build_NEMS.log 2>&1
}

#------------------------------------
# build chgres
#------------------------------------
$Build_chgres && {
echo " .... Building chgres .... "
./build_chgres.sh > $logs_dir/build_chgres.log 2>&1
}

#------------------------------------
# build sfcanl_nsttfchg 
#------------------------------------
$Build_sfcanl_nsttfchg && {
echo " .... Building nst_tf_chg .... "
./build_sfcanl_nsttfchg.sh > $logs_dir/build_sfcanl_nsttfchg.log 2>&1
}

#------------------------------------
# build orog
#------------------------------------
$Build_orog && {
echo " .... Building orog .... "
./build_orog.sh > $logs_dir/build_orog.log 2>&1
}

#------------------------------------
# build cycle 
#------------------------------------
$Build_cycle && {
echo " .... Building cycle .... "
./build_cycle.sh > $logs_dir/build_cycle.log 2>&1
}

#------------------------------------
# build tropcy_NEMS
#------------------------------------
$Build_tropcy && {
echo " .... Building tropcy_NEMS .... "
./build_tropcy_NEMS.sh > $logs_dir/build_tropcy_NEMS.log 2>&1
}

#------------------------------------
# build emcsfc
#------------------------------------
$Build_emcsfc && {
echo " .... Building emcsfc .... "
./build_emcsfc.sh > $logs_dir/build_emcsfc.log 2>&1
}

#------------------------------------
# build fre-nctools
#------------------------------------
$Build_nctools && {
echo " .... Building fre-nctools .... "
./build_fre-nctools.sh > $logs_dir/build_fre-nctools.log 2>&1
}

echo;echo " .... Build system finished .... "

exit 0
