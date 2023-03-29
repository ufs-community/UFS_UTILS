help([[
Load environment to compile UFS_UTILS on Hera using Gnu
]])

cmake_ver=os.getenv("cmake_ver") or "3.16.1"
load(pathJoin("cmake", cmake_ver))

hpss_ver=os.getenv("hpss_ver") or ""
load(pathJoin("hpss", hpss_ver))

prepend_path("MODULEPATH", "/scratch1/NCEPDEV/nems/role.epic/gnu/modulefiles")

gnu_ver="9.2"
load(pathJoin("gnu", gnu_ver))

prepend_path("MODULEPATH", "/scratch1/NCEPDEV/nems/role.epic/hpc-stack/libs/gnu-9.2_ncdf49/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.2.0"
load(pathJoin("hpc", hpc_ver))

load(pathJoin("hpc-gnu", gnu_ver))

mpich_ver=os.getenv("mpich_ver") or "3.3.2"
load(pathJoin("mpich", mpich_ver))

load(pathJoin("hpc-mpich", mpich_ver))

netcdf_ver=os.getenv("netcdf_ver") or "4.9.1"
load(pathJoin("netcdf", netcdf_ver))

esmf_ver=os.getenv("esmf_ver") or "8.4.1"
load(pathJoin("esmf", esmf_ver))

bacio_ver=os.getenv("bacio_ver") or "2.4.1"
load(pathJoin("bacio", bacio_ver))

g2_ver=os.getenv("g2_ver") or "3.4.5"
load(pathJoin("g2", g2_ver))

ip_ver=os.getenv("ip_ver") or "3.3.3"
load(pathJoin("ip", ip_ver))

nemsio_ver=os.getenv("nemsio_ver") or "2.5.4"
load(pathJoin("nemsio", nemsio_ver))

sp_ver=os.getenv("sp_ver") or "2.3.3"
load(pathJoin("sp", sp_ver))

w3nco_ver=os.getenv("w3nco_ver") or "2.4.1"
load(pathJoin("w3nco", w3nco_ver))

sfcio_ver=os.getenv("sfcio_ver") or "1.4.1"
load(pathJoin("sfcio", sfcio_ver))

sigio_ver=os.getenv("sigio_ver") or "2.3.2"
load(pathJoin("sigio", sigio_ver))

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

zlib_ver=os.getenv("zlib_ver") or "1.2.11"
load(pathJoin("zlib", zlib_ver))

png_ver=os.getenv("png_ver") or "1.6.37"
load(pathJoin("libpng", png_ver))

whatis("Description: UFS_UTILS build environment")
