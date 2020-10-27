
# UFS_UTILS

Utilities for the NCEP models. This is part of the
[NCEPLIBS](https://github.com/NOAA-EMC/NCEPLIBS) project.

## Authors

NCEP/EMC developers.

Code manager: George Gayno

## Prerequisites

This package requires the following NCEPLIBS packages:
 - [NCEPLIBS-gfsio](https://github.com/NOAA-EMC/NCEPLIBS-gfsio)
 - [NCEPLIBS-sfcio](https://github.com/NOAA-EMC/NCEPLIBS-sfcio)
 - [NCEPLIBS-w3nco](https://github.com/NOAA-EMC/NCEPLIBS-w3nco)
 - [NCEPLIBS-landsfcutil](https://github.com/NOAA-EMC/NCEPLIBS-landsfcutil)
 - [NCEPLIBS-bacio](https://github.com/NOAA-EMC/NCEPLIBS-bacio)
 - [NCEPLIBS-nemsio](https://github.com/NOAA-EMC/NCEPLIBS-nemsio)
 - [NCEPLIBS-nemsiogfs](https://github.com/NOAA-EMC/NCEPLIBS-nemsiogfs)
 - [NCEPLIBS-sigio](https://github.com/NOAA-EMC/NCEPLIBS-sigio)
 - [NCEPLIBS-sp](https://github.com/NOAA-EMC/NCEPLIBS-sp)
 - [NCEPLIBS-ip](https://github.com/NOAA-EMC/NCEPLIBS-ip)
 - [NCEPLIBS-w3emc](https://github.com/NOAA-EMC/NCEPLIBS-w3emc)
 - [NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2)
 - [NCEPLIBS-wgrib2](https://github.com/NOAA-EMC/NCEPLIBS-wgrib2)

This package also requires:

 - [netcdf-c Library](https://github.com/Unidata/netcdf-c)
 - [netcdf-fortran Library](https://github.com/Unidata/netcdf-fortran)
 - [ESMF](https://github.com/esmf-org/esmf)

## Installing

```
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=/path/to/install ..
make -j2
make install
```

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

