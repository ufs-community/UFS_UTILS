@mainpage

# UFS_UTILS

Utilities for the NCEP models. This is part of the
[NCEPLIBS](https://github.com/NOAA-EMC/NCEPLIBS) project.

The UFS_UTILS code can be found here:
https://github.com/NOAA-EMC/UFS_UTILS.

## The Utilities

- chgres_cube - Creates cold start initial conditions for FV3 model
  runs.

- emcsfc_ice_blend - Blends National Ice Center sea ice cover and EMC
  sea ice concentration data to create a global sea ice analysis used
  to update the GFS once per day.

- emcsfc_snow2mdl - Blends National Ice Center snow cover and Air
  Force snow depth data to create a global depth analysis used to
  update the GFS snow field once per day.

- fre-nctools - Tools to remap data; and to create the geo-reference
  fields (latitude, longitude, etc.) for an FV3 grid.

- fvcom_tools - - Replaces lake surface and lake ice temperature along
  with aerial ice concentration generated from the Great Lakes
  Operational Forecast System (GLOFS) in an FV3 surface restart
  file. [fvcom documentation](sorc/fvcom_tools.fd/fvcom_readme.md)
 
- global_chgres - Creates cold start initial conditions for FV3 model
  runs. Deprecated by the chgres_cube utility. [global_chgres
  documentation](sorc/global_chgres.fd/global_chgres_users_guide.md)

- global_cycle - Updates the GFS surface conditions using external
  snow and sea ice analyses. Updates monthly climatological fields
  such as plant greenness fraction and albedo. Runs as part of the GFS
  and GDAS cycles.

- grid_tools - Utilities to filter topography, to create regional
  extended Schmidt gnomonic grids, and to compute the equivalent
  global resolution of a regional grid.

- nst_tf_chg

- orog_mask_tools

- sfc_climo_gen

- vcoord_gen

