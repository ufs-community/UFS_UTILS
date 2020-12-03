#!/bin/bash

#BSUB -W 0:01
#BSUB -o log
#BSUB -e log
#BSUB -J vcoord
#BSUB -q debug
#BSUB -R "rusage[mem=100]"
#BSUB -P GFS-DEV

#-------------------------------------------------------------------------------
#
# Generate a hybrid coordinate interface profile on WCOSS-Cray.
#
# Build the repository using the ./build_all.sh script before running.
#
# Output 'ak' and 'bk' values are placed in $outfile.
#
#-------------------------------------------------------------------------------

set -x

source ../../sorc/machine-setup.sh > /dev/null 2>&1
module use ../../modulefiles
module load build.$target.intel
module list

outfile="./global_hyblev.txt"

levs=128               # integer number of levels
lupp=88                # integer number of levels below pupp
pbot=100000.0          # real nominal surface pressure (Pa)
psig=99500.0           # real nominal pressure where coordinate changes
                       # from pure sigma (Pa)
ppre=7000.0            # real nominal pressure where coordinate changes
                       # to pure pressure (Pa)
pupp=7000.0            # real nominal pressure where coordinate changes
                       # to upper atmospheric profile (Pa)
ptop=0.0               # real pressure at top (Pa)
dpbot=240.0            # real coordinate thickness at bottom (Pa)
dpsig=1200.0           # real thickness of zone within which coordinate changes
                       # to pure sigma (Pa)
dppre=18000.0          # real thickness of zone within which coordinate changes
                       # to pure pressure (Pa)
dpupp=550.0            # real coordinate thickness at pupp (Pa)
dptop=1.0              # real coordinate thickness at top (Pa)

rm -f $outfile

echo $levs $lupp $pbot $psig $ppre $pupp $ptop $dpbot $dpsig $dppre $dpupp $dptop | $PWD/../../exec/vcoord_gen > $outfile

exit
