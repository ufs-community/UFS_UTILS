
# regrid_sfc

# Introduction

The program regridStates.F90 is used to regrid masked sfc fields
using masked bi-linear interpolaton. It is designed for 
regridding of increments and states, and can regrid from/to 
the FV3 grid that the GFS model runs on, and the Gaussian grid 
used for GFS model output (i.e., history files).

Clara Draper, Aug 17, 2024.

* Uses bi-linear interpolation, with masking of the input and output. 
* Current masking options are all soil (land, no glaciers) and snow-free soil (input increment only)
* Does not guarantee all output grid cells will have mapped values, by design. CANNOT BY USED FOR RESTARTS (could be changed though).
* Output grid cells that do not have a value after the interpolation step are filled with the nearest neighbouring output grid cell. If there are no mapped neighbours within two layers of surrounding neighbours, the output grid cell will remain unmapped. In this way, islands do not get filled with distant values. 
* Reads in / writes out fv3 files and Gaussian grid (gsi output increment) files.
* Handling Gaussian files requires a scrip file with the grid details. This can be generated using the ufs_utils weigh_gen program.
* has separate namelist for input and output. Required options depend on the grid type and whether it's in/out. Look in readin_setup routine to check what is needed. 

FEATURES TO ADD:
* Can only process 2D variables (no vertical dimension; TO-DO fix this) 
* Re-think parallelization so can do multiple ensemble members at once
* Add code to check required options are available for requested set-up grid

This document is part of the <a href="../index.html">UFS_UTILS
documentation</a>.

The regridStates program is part of the
[UFS_UTILS](https://github.com/ufs-community/UFS_UTILS) project.
