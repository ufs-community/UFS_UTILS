#!/bin/bash

set -x

levs=128
lupp=88
pbot=100000.0
psig=99500.0
ppre=7000.0
pupp=7000.0
ptop=0.0
dpbot=240.0
dpsig=1200.0
dppre=18000.0
dpupp=550.0
dptop=1.0

#echo 128 88 100000. 99500. 7000. 7000. 0. 240. 1200. 18000. 550. 1.0 | ../../exec/wcoss_ttakbkgen > global_hyblev.txt
echo $levs $lupp $pbot $psig $ppre $pupp $ptop $dpbot $dpsig $dppre $dpupp $dptop | ../../exec/wcoss_ttakbkgen > global_hyblev.txt
