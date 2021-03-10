@mainpage

# UFS_UTILS

Utilities for the NCEP models. This is part of the
[NCEPLIBS](https://github.com/NOAA-EMC/NCEPLIBS) project.

The UFS_UTILS code can be found here:
https://github.com/NOAA-EMC/UFS_UTILS.

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


