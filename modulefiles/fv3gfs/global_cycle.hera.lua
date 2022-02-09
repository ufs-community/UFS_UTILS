help([[
Load environment to compile UFS_UTILS global_cycle on Hera
]])

prepend_path("MODULEPATH", "/scratch2/NCEPDEV/nwprod/hpc-stack/libs/hpc-stack/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.1.0"
load(pathJoin("hpc", hpc_ver))

hpc_intel_ver=os.getenv("hpc_intel_ver") or "18.0.5.274"
load(pathJoin("hpc-intel", hpc_intel_ver))

hpc_intel_ver=os.getenv("hpc_impi_ver") or "2018.0.4"
load(pathJoin("hpc-impi", hpc_impi_ver))

ip_ver=os.getenv("ip_ver") or "3.3.3"
load(pathJoin("ip", ip_ver))

sp_ver=os.getenv("sp_ver") or "2.3.3"
load(pathJoin("sp", sp_ver))

w3nco_ver=os.getenv("w3nco_ver") or "2.4.1"
load(pathJoin("w3nco", w3nco_ver))

bacio_ver=os.getenv("bacio_ver") or "2.4.1"
load(pathJoin("bacio", bacio_ver))

prepend_path("MODULEPATH", "/scratch2/NCEPDEV/nwprod/NCEPLIBS/modulefiles")
netcdf_ver=os.getenv("netcdf_ver") or "4.7.4"
load(pathJoin("netcdf_parallel", netcdf_ver))

setenv("FCMP","mpiifort")

whatis("Description: UFS_UTILS global_cycle build environment")
