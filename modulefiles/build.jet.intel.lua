help([[
Load environment to compile UFS_UTILS on Jet
]])

cmake_ver=os.getenv("cmake_ver") or "3.16.1"
load(pathJoin("cmake", cmake_ver))

hpss_ver=os.getenv("hpss_ver") or ""
load(pathJoin("hpss", hpss_ver))

prepend_path("MODULEPATH", "/lfs4/HFIP/hfv3gfs/role.epic/hpc-stack/libs/intel-2022.1.2/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.2.0"
load(pathJoin("hpc", hpc_ver))

hpc_intel_ver=os.getenv("hpc_intel_ver") or "2022.1.2"
load(pathJoin("hpc-intel", hpc_intel_ver))

impi_ver=os.getenv("impi_ver") or "2022.1.2"
load(pathJoin("hpc-impi", impi_ver))

hdf5_ver=os.getenv("hdf5_ver") or "1.10.6"
load(pathJoin("hdf5", hdf5_ver))

netcdf_ver=os.getenv("netcdf_ver") or "4.7.4"
load(pathJoin("netcdf", netcdf_ver))

nccmp_ver=os.getenv("nccmp_ver") or "1.8.9.0"
load(pathJoin("nccmp", nccmp_ver))

esmf_ver=os.getenv("esmf_ver") or "8.4.0b08"
load(pathJoin("esmf", esmf_ver))

w3nco_ver=os.getenv("w3nco_ver") or "2.4.1"
load(pathJoin("w3nco", w3nco_ver))

sp_ver=os.getenv("sp_ver") or "2.3.3"
load(pathJoin("sp", sp_ver))

ip_ver=os.getenv("ip_ver") or "3.3.3"
load(pathJoin("ip", ip_ver))

bacio_ver=os.getenv("bacio_ver") or "2.4.1"
load(pathJoin("bacio", bacio_ver))

sigio_ver=os.getenv("sigio_ver") or "2.3.2"
load(pathJoin("sigio", sigio_ver))

sfcio_ver=os.getenv("sfcio_ver") or "1.4.1"
load(pathJoin("sfcio", sfcio_ver))

nemsio_ver=os.getenv("nemsio_ver") or "2.5.4"
load(pathJoin("nemsio", nemsio_ver))

g2_ver=os.getenv("g2_ver") or "3.4.5"
load(pathJoin("g2", g2_ver))

prod_util_ver=os.getenv("prod_util_ver") or "1.2.2"
load(pathJoin("prod_util", prod_util_ver))

nco_ver=os.getenv("nco_ver") or "4.9.3"
load(pathJoin("nco", nco_ver))

whatis("Description: UFS_UTILS build environment")
