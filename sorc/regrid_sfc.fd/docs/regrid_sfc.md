
# regrid_sfc

# Introduction

The program regridStates.F90 is used to regrid masked sfc fields
using masked bi-linear interpolaton. It is designed for 
regridding of increments and states, and can regrid from/to 
the FV3 grid that the GFS model runs on, and the Gaussian grid 
used for GFS model output (i.e., history files).

This document is part of the <a href="../index.html">UFS_UTILS
documentation</a>.

The regridStates program is part of the
[UFS_UTILS](https://github.com/ufs-community/UFS_UTILS) project.
