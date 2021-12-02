help([[
Load environment to compile UFS_UTILS global_cycle on Jet
]])

intel_ver=os.getenv("intel_ver") or "18.0.5.274"
load(pathJoin("intel", intel_ver))

mpich_intel_ver=os.getenv("mpich_intel_ver") or "2018.4.274"
load(pathJoin("impi", mpich_intel_ver))

w3nco_ver=os.getenv("w3nco_ver") or "v2.0.6"
load(pathJoin("w3nco", w3nco_ver))

sp_ver=os.getenv("sp_ver") or "v2.0.2"
load(pathJoin("sp", sp_ver))

bacio_ver=os.getenv("bacio_ver") or "v2.0.2"
load(pathJoin("bacio", bacio_ver))

ip_ver=os.getenv("ip_ver") or "v3.0.0"
load(pathJoin("ip", ip_ver))

szip_ver=os.getenv("szip_ver") or ""
load(pathJoin("szip", szip_ver))

hdf5_ver=os.getenv("hdf5_ver") or "1.8.9"
load(pathJoin("hdf5", hdf5_ver))

netcdf_ver=os.getenv("netcdf_ver") or "4.2.1.1"
load(pathJoin("netcdf", netcdf_ver))
setenv("NETCDF_INCLUDE","-I${NETCDF}/include")
setenv("NETCDF_LDFLAGS_F","-L${NETCDF}/lib -lnetcdf -lnetcdff -L${HDF5}/lib -lhdf5 -lhdf5_fortran")

setenv("FCMP","mpiifort")

whatis("Description: UFS_UTILS global_cycle build environment")
