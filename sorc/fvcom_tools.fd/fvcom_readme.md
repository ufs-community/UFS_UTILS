@brief replaces lake surface and lake ice temperature @anchor fvcom_readme

**fvcom_to_FV3.exe**

**Introduction:**
 This code replaces lake surface and lake ice temperature along
 with aerial ice concentration generated from Great Lakes 
 Operational Forecast System (GLOFS), an FVCOM-based model, into 
 sfc_data.nc.
 **NOTE** that the variables in the input files must reside on 
 the same grid. This means data from FVCOM must be horizontally 
 interpolated to the FV3 grid. This routine will also force a 
 minimum ice concentration of 15%. If ice concentration is less 
 than 15% in FVCOM, it will be set to 0% to avoid FV3 from 
 changing values less than 15% to 15% and generating unrealistic 
 lake ice temperatures.

**Library Dependencies:**
 Installation depends on the netCDF library and cmake.

**Running:**
 This routine will take two variables from the command line:
 1. Name of FV3 sfc data file (e.g. sfc_data.tile7.halo0.nc)
   which is generated from chgres_cube.exe.
 2. Name of FVCOM data file in netcdf format (e.g. fvcom.nc)

 To run the script, use the following example, modifying file
 names as needed:
   ./fvcom_to_FV3 sfc_data.tile7.halo0.nc fvcom.nc
 Output will be to the sfc data file and include lake surface 
 and lake ice temperature, and lake ice concentration from FVCOM.


This routine is *strongly* based upon Eric James' (ESRL/GSL) work
 to update HRRR/WRF Great Lakes' temperature data with FVCOM.
 It also relies heavily on Ming Hu's (ESRL/GSL) ncio module.

**For more information, please contact:**
 David Wright
 University of Michigan and GLERL
 dmwright@umich.edu
