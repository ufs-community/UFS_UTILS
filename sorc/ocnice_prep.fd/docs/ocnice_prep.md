# ocnice_prep

# Introduction

The ocnice_prep program will downscale a MOM6 or CICE6 1/4deg tripole restart file to a desired lower resolution tripole warm-start file using ESMF regridding.

This document is part of the <a href="../index.html">UFS_UTILS
documentation</a>.

The ocnice_prep program is part of the
[UFS_UTILS](https://github.com/ufs-community/UFS_UTILS) project.

## Creating a warm-start file for MOM6 or CICE6 from a restart file at higher resolution

A warmstart file can be created from an existing restart file (at higher resolution) by using ESMF Regridding to map fields
to a destination grid. For MOM6, the procedure produces a warm-start file (composed of T,S,U,V and Interface Height), where the
utility creates the interface heights (``eta``) using the sea surface height (``sfc``) and the interface thickness (``h``) in
the MOM6 restart. A full restart file is not generated because MOM6 has resolution dependent parameterizations, so that restart
files for different resolutions might contain different fields. For CICE6, a file consisting of all the restart fields is produced.
While in the form of a "true" restart, it is also more appropriately considered to be a warm-start file.

For MOM6, the following MOM_input settings can then be used with the warm-start file (e.g. ``ocean.mx100.nc``):

```
! === module MOM_state_initialization ===
INIT_LAYERS_FROM_Z_FILE = False  !   [Boolean] default = False
                                ! If true, initialize the layer thicknesses, temperatures, and salinities from a
                                ! Z-space file on a latitude-longitude grid.
! WARMSTARTS
THICKNESS_CONFIG = "file"       ! default = "uniform"
                                ! A string that determines how the initial layer thicknesses are specified for a
                                ! new run:
                                !     file - read interface heights from the file specified
                                !       by (THICKNESS_FILE).
THICKNESS_FILE = "ocean.mx100.nc" !
                                ! The name of the thickness file.
ADJUST_THICKNESS = True         !   [Boolean] default = False
                                ! If true, all mass below the bottom removed if the topography is shallower than
                                ! the thickness input file would indicate.
INTERFACE_IC_VAR = "eta"        ! default = "eta"
                                ! The variable name for initial conditions for interface heights relative to
                                ! mean sea level, positive upward unless otherwise rescaled.
TS_CONFIG = "file"              !
                                ! A string that determines how the initial temperatures and salinities are
                                ! specified for a new run:
                                !     file - read velocities from the file specified
TS_FILE = "ocean.mx100.nc"      !
                                ! The initial condition file for temperature.
TEMP_IC_VAR = "Temp"            ! default = "PTEMP"
                                ! The initial condition variable for potential temperature.
SALT_IC_VAR = "Salt"            ! default = "SALT"
                                ! The initial condition variable for salinity.
SALT_FILE = "ocean.mx100.nc"    ! default = "ocean.mx100.nc"
                                ! The initial condition file for salinity.
VELOCITY_CONFIG = "file"        ! default = "zero"
                                ! A string that determines how the initial velocities are specified for a new
                                ! run:
                                !     file - read velocities from the file specified
                                !       by (VELOCITY_FILE).
VELOCITY_FILE = "ocean.mx100.nc" !
                                ! The name of the velocity initial condition file.
U_IC_VAR = "u"                  ! default = "u"
                                ! The initial condition variable for zonal velocity in VELOCITY_FILE.
V_IC_VAR = "v"                  ! default = "v"
```


For CICE6, the warm-start file (e.g. ``ice.mx100.nc``) can be used directly in the ``ice_in`` namelist:

```
    runtype             = 'initial'
    runid               = 'unknown'
    ice_ic              = 'ice.mx100.nc'
```

## Remapping procedure

The remapping is done using an ESMF remapping via a RouteHandle, where the both the source and destination mask values are ``0``.
This ensures that only water points are mapped between the grids. Extrapolation is used to map destination points which are unmapped
because the source grid was not a water point. Since for the ocean, the required mapping varies with depth, the utility makes
use of dynamic masking if required. This allows the same RouteHandle to be used for all depths, with masking at each depth done
_on the fly_.

Mapping for MOM6 is done using ``bilinear`` mapping, whereas for CICE6 it is done using ``nearest-source-to-destination`` mapping.
For CICE6, this ensures that the thermodynamic fields in the CICE6 restart, which are for each vertical layer (typically 7)
and each thickness category (typically 5), are consistent when mapped.

Fields in the restart files for both MOM6 and CICE6 are defined at locations on the Arakawa C-grid. Similar to the ``cpld_gridgen``
utility, these are referred to here as the ``Ct``, ``Cu``, ``Cv`` and ``Bu`` locations. Vectors in MOM6 are natively located at
the north and east faces (``Cv`` and ``Cu``) while CICE vectors are located at ``Bu`` locations. **NOTE:** The CICE-C grid is
not currently supported.

The general procedure, for either MOM6 or CICE6 is as follows:

1. Fields are retrieved from the source file and placed on the ``Ct`` grid location. For vector fields, the vectors are rotated
from their orientation along model dimensions (IJ) to eastward-northward (EN) directions.
2. Fields are packed into arrays by mapping type (bilinear or conservative) and dimensionality.
3. Fields are remapped using ESMF. This remapping may use dynamic Masking.
4. Vector fields are re-rotated back to model index direction (EN->IJ) and remapped back to their native stagger locations.
5. Fields are written into a warm-start file.

## Required files

The following files are required.

- The ESMF mesh file for the source and destination grids. These are available as products of the ``cpld_gridgen`` utility.
- The remapping weights for mapping fields to and from ``Ct`` grid locations on both the souce and destination grids. These are also available as products of the ``cpld_gridgen`` utility.
- A text file (csv) listing the fields to be remapped for either the ocean or ice restart file.
- A namelist file defining the type of regridding desired (``ocean`` or ``ice``) and the locations of the needed weights.

## Pre-generation of MOM6 restart file

Because the required fields for MOM6 are located in two separate restart files, a single file containing all the necessary fields
must be generated for the ``ocean`` case. This is done using NetCDF operators (NCO) with the following commands, assuming that the
two required MOM6 restart files are available locally and are named ``MOM.res.nc`` and ``MOM.res_1.nc``.

```
ncks -v Temp,Salt,h,u MOM.res.nc ocean.nc
ncks -v v,sfc -A MOM.res_1.nc ocean.nc
```
