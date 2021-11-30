help([[
Load environment to compile UFS_UTILS global_cycle on WCOSS-Dell
]])

intel_ver=os.getenv("intel_ver") or "18.0.1.163"
load(pathJoin("ips", intel_ver))

cray_mpich_ver=os.getenv("cray_mpich_ver") or "18.0.1"
load(pathJoin("impi", cray_mpich_ver))

netcdf_ver=os.getenv("netcdf_ver") or "4.7.4"
load(pathJoin("NetCDF-parallel", netcdf_ver))

w3nco_ver=os.getenv("w3nco_ver") or "2.0.6"
load(pathJoin("w3nco", w3nco_ver))

sp_ver=os.getenv("sp_ver") or "2.0.2"
load(pathJoin("sp", sp_ver))

ip_ver=os.getenv("ip_ver") or "3.0.1"
load(pathJoin("ip", ip_ver))

bacio_ver=os.getenv("bacio_ver") or "2.0.2"
load(pathJoin("bacio", bacio_ver))

setenv("FCMP","mpif90")

whatis("Description: UFS_UTILS global_cycle build environment")
