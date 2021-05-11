
# noah

# Introduction

This directory contains copies of the NOAH LSM code needed to make any land state adjusments made neceessary by the DA updates.
For the NOAH model, in the current code routines not labelled public/private, will shorty be replaced, and are very unlikely to change. For now, have made a second copy of the relevant routines.

The files were checked out from https://github.com/NCAR/ccpp-physics -b master, commit 08b72bc1c23c48a823626d81f8e0a398685a35a3 (dated May, 2021).

Files used from the checkout:
namelist_soilveg.f
set_soilveg.f      
sflx.F
machine.F
physcons.F90

sflx_snippet.f is a snippet of the file sflx.F, with minor changes to compile, and to enable the frh2o routine to be called externally. Copied values from physcons and machine into sflx_snippet. 

This document is part of the <a href="../index.html">UFS_UTILS documentation.</a> 

The NOAH routines are included in the <a href="../lsm_routines/index.html">lsm_routines directory.</a>

The NOAH library created here is used in the  [NCEPLIBS
UFS_UTILS](https://github.com/NOAA-EMC/UFS_UTILS) project.
