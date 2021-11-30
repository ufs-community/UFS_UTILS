help([[
Load environment to compile UFS_UTILS global_cycle on Orion
]])

intel_ver=os.getenv("intel_ver") or "2020"
load(pathJoin("intel", intel_ver))

cray_mpich_ver=os.getenv("cray_mpich_ver") or "2020"
load(pathJoin("impi", cray_mpich_ver))

prepend_path("MODULEPATH", "/apps/contrib/NCEPLIBS/orion/modulefiles")

w3nco_ver=os.getenv("w3nco_ver") or "2.1.0"
load(pathJoin("w3nco", w3nco_ver))

ip_ver=os.getenv("ip_ver") or "3.1.0"
load(pathJoin("ip", ip_ver))

sp_ver=os.getenv("sp_ver") or "2.0.3"
load(pathJoin("sp", sp_ver))

bacio_ver=os.getenv("bacio_ver") or "2.2.0"
load(pathJoin("bacio", bacio_ver))

prepend_path("MODULEPATH", "/apps/contrib/NCEPLIBS/lib/modulefiles")

netcdf_ver=os.getenv("netcdf_ver") or "4.7.4.release"
load(pathJoin("netcdfp", netcdf_ver))
setenv("NETCDF_INCLUDE","${NETCDF_FFLAGS}")
setenv("NETCDF_LDFLAGS_F","-L${NETCDF}/lib -lnetcdf -lnetcdff -L${HDF5}/lib -lhdf5")

setenv("FCMP","mpiifort")

whatis("Description: UFS_UTILS global_cycle build environment")
