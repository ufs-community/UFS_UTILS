#!/bin/bash

#-------------------------------------------------------------------------------
#
# Generate a hybrid coordinate interface profile.
# On WCOSS2, use 'run.wcoss2.sh'.
#
# Build the repository using the ./build_all.sh script before running.
#
# Output 'ak' and 'bk' values are placed in $outfile.
#
#-------------------------------------------------------------------------------

set -x

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

echo $levs $lupp $pbot $psig $ppre $pupp $ptop $dpbot $dpsig $dppre $dpupp $dptop | ../../exec/vcoord_gen > $outfile

exit
