# fvcom_tools

# Introduction

This code replaces lake surface and lake ice temperature along
with aerial ice concentration generated from Great Lakes 
Operational Forecast System (GLOFS), an FVCOM-based model, into 
sfc_data.nc.

This document is part of the <a href="../index.html">UFS_UTILS
documentation</a>.

The fvcom_tools program is part of the
[UFS_UTILS](https://github.com/ufs-community/UFS_UTILS) project.

## NOTE

The variables in the input files must reside on 
the same grid. This means data from FVCOM must be horizontally 
interpolated to the FV3 grid. This routine will also force a 
minimum ice concentration of 15%. If ice concentration is less 
than 15% in FVCOM, it will be set to 0% to avoid FV3 from 
changing values less than 15% to 15% and generating unrealistic 
lake ice temperatures.

## Library Dependencies:

Installation depends on the netCDF library and cmake.

## Running

This routine will take four variables from the command line:
1. Name of FV3 sfc data file (e.g. sfc_data.tile7.halo0.nc)
   which is generated from chgres_cube.exe.
2. Name of FVCOM data file in netcdf format (e.g. fvcom.nc)
3. "warm" or "cold" start. "warm" start will read in 
    sfc_data.nc files generated from a restart of UFS-SRW.
    "cold" start will read in sfc_data.nc files generated 
    from chgres_cube. 
4. String of time slice to use in the fvcom.nc file. This string
    should match exactly what is in the Times variable of the .nc file.
To run the script, use the following example, modifying file
names as needed:
   ./fvcom_to_FV3 sfc_data.tile7.halo0.nc fvcom.nc cold \
     2020-01-31T18:00:00.000000
Output will be to the sfc data file and include lake surface 
and lake ice temperature, and lake ice concentration from the 
first time in the FVCOM file.

This routine is *strongly* based upon Eric James' (ESRL/GSL) work
to update HRRR/WRF Great Lakes' temperature data with FVCOM.
It also relies heavily on Ming Hu's (ESRL/GSL) ncio module.

## For more information, please contact:

David Wright, University of Michigan and GLERL: dmwright@umich.edu
