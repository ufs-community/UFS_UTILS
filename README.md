
# UFS_UTILS

Utilities for the NCEP models. This is part of the
[NCEPLIBS](https://github.com/NOAA-EMC/NCEPLIBS) project.

Documentation for chgres_cube and other utilities can be found at
https://noaa-emcufs-utils.readthedocs.io/en/latest/.

Complete documentation can be found at
https://ufs-community.github.io/UFS_UTILS/.

## Authors

Utility | Programmer(s)
--------|----------
chgres_cube | George Gayno, Jeff Beck, Larissa Reames
cpld_gridgen | Denise Worthen
emcsfc_ice_blend | George Gayno
emcsfc_snow2mdl | George Gayno
fre-nctools | GFDL progammer
fvcom_tools | David Wright, University of Michigan, Ming Hu, GSD/AMB
gblevents | Hang Lei
gdas_init | George Gayno
global_cycle | George Gayno, Shrinivas Moorthi, Mark Iredell, Xu Li, Hang Lei
grid_tools | R. J. Purser (regional_esg_grid), Tom Black/Ben Blake (shave.fd), Gerard Ketefian (global_equiv_resol), Tsukasa Fujita, JMA (pmat2), GFDL programmer (topo filtering code).
orog_mask_tools | Ning Wang, Jordan Alpert, Shan Sun and Ning Wang
sfc_climo_gen | George Gayno
vcoord_gen | Fanglin Yang, Mark Iredell

UFS_UTILS Code managers: George Gayno, Kyle Gerheiser, Jeff Beck, Larissa Reames

## Prerequisites

This package uses the [hpc-stack](https://github.com/NOAA-EMC/hpc-stack) for the following NCEPLIBS packages:
 - [NCEPLIBS-sfcio](https://github.com/NOAA-EMC/NCEPLIBS-sfcio)
 - [NCEPLIBS-w3nco](https://github.com/NOAA-EMC/NCEPLIBS-w3nco)
 - [NCEPLIBS-bacio](https://github.com/NOAA-EMC/NCEPLIBS-bacio)
 - [NCEPLIBS-nemsio](https://github.com/NOAA-EMC/NCEPLIBS-nemsio)
 - [NCEPLIBS-sigio](https://github.com/NOAA-EMC/NCEPLIBS-sigio)
 - [NCEPLIBS-sp](https://github.com/NOAA-EMC/NCEPLIBS-sp)
 - [NCEPLIBS-ip](https://github.com/NOAA-EMC/NCEPLIBS-ip)
 - [NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2)

And for the following third party libraries:

 - [netcdf-c Library](https://github.com/Unidata/netcdf-c)
 - [netcdf-fortran Library](https://github.com/Unidata/netcdf-fortran)
 - [ESMF](https://github.com/esmf-org/esmf)
 - [Jasper](https://github.com/jasper-software/jasper)
 - [Zlib](www.zlib.net)
 - [HDF5](https://www.hdfgroup.org/solutions/hdf5/)
 - [PNG](http://www.libpng.org/pub/png/)

It also uses the following repositories:

 - [NCAR common community physics package](https://github.com/NCAR/ccpp-physics)

## Installing

On Orion, Jet, Hera and WCOSS2, invoke the build script:

```
./build_all.sh
```

Otherwise, do:

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
make -j2
make install
```

## Contents

The UFS_UTILS package contains the following utilities (under the sorc
directory):
- chgres_cube
- cpld_gridgen
- emcsfc_ice_blend
- emcsfc_snow2mdl
- fre-nctools
- fvcom_tools
- gblevents
- global_cycle
- grid_tools
- orog_mask_tools
- sfc_climo_gen
- vcoord_gen

The reg_tests directory contains the consistency test code.

The fix directory is where we set links to directories containing
large, static data files used by UFS_UTILS programs.

The tests directory contains unit tests.

The ush directory contains scripts to run UFS_UTILS programs.  Most
are called from driver scripts.

The util directory contains utility scripts to create coldstart
initial conditions for GFS parallels, and to run the vertical
coordinate generator.

The parm directory contains variable mapping parameter tables used by
the chgres_cube program.

The driver_scripts directory contains high-level driver scripts to
create a model grid on officially supported HPC platforms.

The modulefiles directory contains modules loaded when building
UFS_UTILS on supported HPC platforms.  They are also loaded at runtime
by utility and consistency test scripts.

The docs directory contains the control file for the doxygen
documentation build, as well as some markdown files which are part of
the documentation. It also contains (in the source subdirectory) the
ReadTheDocs documentation files.

The cmake directory contains CMake package find utilities, and utilities to
run units tests on some supported HPC platforms.

## References

Gayno G., Beck J., Carson L., [Pre-Processing:
chgres_cube](./docs/20201105-0945a-pre-processing-chgres-cube-gayno-final.pdf),
UFS MRW App Training, 5 November 2020.

## Disclaimer

The United States Department of Commerce (DOC) GitHub project code is
provided on an "as is" basis and the user assumes responsibility for
its use. DOC has relinquished control of the information and no longer
has responsibility to protect the integrity, confidentiality, or
availability of the information. Any claims against the Department of
Commerce stemming from the use of its GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial
products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of
Commerce. The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.

