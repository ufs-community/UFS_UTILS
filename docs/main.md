@mainpage

# UFS_UTILS

Utilities for the NCEP Unified Forecast System.

The UFS_UTILS code can be found here:
https://github.com/ufs-community/UFS_UTILS.

## Documentation for Previous Versions of UFS_UTILS

* [UFS_UTILS Version 1.11.0](ver-1.11.0/index.html)
* [UFS_UTILS Version 1.10.0](ver-1.10.0/index.html)
* [UFS_UTILS Version 1.9.0](ver-1.9.0/index.html)
* [UFS_UTILS Version 1.8.0](ver-1.8.0/index.html)
* [UFS_UTILS Version 1.7.0](ver-1.7.0/index.html)
* [UFS_UTILS Version 1.6.0](ver-1.6.0/index.html)
* [UFS_UTILS Version 1.5.0](ver-1.5.0/index.html)
* [UFS_UTILS Version 1.4.0](ver-1.4.0/index.html)
* [UFS_UTILS Version 1.3.0](ver-1.3.0/index.html)

## The Utilities

- <a href="chgres_cube/index.html">chgres_cube</a> - Creates cold
  start initial conditions for FV3 model runs.

- <a href="emcsfc_ice_blend/index.html">emcsfc_ice_blend</a> - Blends
  National Ice Center sea ice cover and EMC sea ice concentration data
  to create a global sea ice analysis used to update the GFS once per
  day.

- <a href="emcsfc_snow2mdl/index.html">emcsfc_snow2mdl</a> - Blends
  National Ice Center snow cover and Air Force snow depth data to
  create a global depth analysis used to update the GFS snow field
  once per day.

- fre-nctools - Tools to remap data; and to create the geo-reference
  fields (latitude, longitude, etc.) for an FV3 grid.

- <a href="fvcom_tools/index.html">fvcom_tools</a> - Replaces lake
  surface and lake ice temperature along with aerial ice concentration
  generated from the Great Lakes Operational Forecast System (GLOFS)
  in an FV3 surface restart file.

- <a href="gblevents/index.html">gblevents</a> -
  Prepares observational prepbufr reports for subsequent quality
  control and analysis programs.
 
- <a href="global_cycle/index.html">global_cycle</a> -
  Updates the GFS surface conditions using external snow and sea ice
  analyses. Updates monthly climatological fields such as plant
  greenness fraction and albedo. Runs as part of the GFS and GDAS
  cycles.

- <a href="grid_tools/index.html">grid_tools</a> -
  Utilities to filter topography, to create regional extended Schmidt
  gnomonic grids, and to compute the equivalent global resolution of a
  regional grid.

- <a href="orog_mask_tools/index.html">orog_mask_tools</a> - Utilities
  to create land mask, terrain and gravity wave drag fields; set lake
  fraction and depth; creates an inland land mask.

- <a href="sfc_climo_gen/index.html">sfc_climo_gen</a> - Creates
  surface climatological fields, such as vegetation type and albedo,
  for an FV3 grid.

- <a href="vcoord_gen/index.html">vcoord_gen</a> - Generates hybrid
  coordinate parameters from fields such as surface pressure, model
  top and the number of vertical levels. Outputs the 'ak' and 'bk'
  parameters used by the forecast model to define the hybrid levels.

- <a href="lsm_routines/index.html">lsm_routines</a> - Land surface 
  model-specific routines that are utilised elsewhere within UFS_UTILS.
  Currently, contains the routines required by global_cycle to 
  perform data assimilation updates to land model states

- <a href="cpld_gridgen/index.html">cpld_gridgen</a> - Utility to 
  create the Fix and IC files for the S2SW and S2S applications 

- <a href="weight_gen/index.html">weight_gen</a> - Utility to 
  create gaussian grid ESMF 'scrip' files for use in creating
  ESMF interpolation weight files.
