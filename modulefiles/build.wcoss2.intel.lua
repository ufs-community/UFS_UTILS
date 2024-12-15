help([[
Load environment to compile UFS_UTILS on WCOSS2
]])

cmake_ver=os.getenv("cmake_ver") or "3.20.2"
load(pathJoin("cmake", cmake_ver))

PrgEnv_intel_ver=os.getenv("PrgEnv_intel_ver") or "8.3.3"
load(pathJoin("PrgEnv-intel", PrgEnv_intel_ver))

craype_ver=os.getenv("craype_ver") or "2.7.17"
load(pathJoin("craype", craype_ver))

intel_ver=os.getenv("intel_ver") or "19.1.3.304"
load(pathJoin("intel", intel_ver))

cray_mpich_ver=os.getenv("cray_mpich_ver") or "8.1.12"
load(pathJoin("cray-mpich", cray_mpich_ver))

libjpeg_ver=os.getenv("libjpeg_ver") or "9c"
load(pathJoin("libjpeg", libjpeg_ver))

zlib_ver=os.getenv("zlib_ver") or "1.2.11"
load(pathJoin("zlib", zlib_ver))

libpng_ver=os.getenv("libpng_ver") or "1.6.37"
load(pathJoin("libpng", libpng_ver))

hdf5_ver=os.getenv("hdf5_ver") or "1.14.0"
load(pathJoin("hdf5-C", hdf5_ver))

netcdf_ver=os.getenv("netcdf_ver") or "4.9.2"
load(pathJoin("netcdf-C", netcdf_ver))

bacio_ver=os.getenv("bacio_ver") or "2.4.1"
load(pathJoin("bacio", bacio_ver))

-- Uncomment when CHGRES_ALL is ON.
--sfcio_ver=os.getenv("sfcio_ver") or "1.4.1"
--load(pathJoin("sfcio", sfcio_ver))

w3emc_ver=os.getenv("w3emc_ver") or "2.9.2"
load(pathJoin("w3emc", w3emc_ver))

nemsio_ver=os.getenv("nemsio_ver") or "2.5.4"
load(pathJoin("nemsio", nemsio_ver))

sigio_ver=os.getenv("sigio_ver") or "2.3.2"
load(pathJoin("sigio", sigio_ver))

sp_ver=os.getenv("sp_ver") or "2.4.0"
load(pathJoin("sp", sp_ver))

ip_ver=os.getenv("ip_ver") or "4.0.0"
load(pathJoin("ip", ip_ver))

g2_ver=os.getenv("g2_ver") or "3.4.5"
load(pathJoin("g2", g2_ver))

-- Needed for mpiexec command.
cray_pals_ver=os.getenv("cray_pals_ver") or "1.2.2"
load(pathJoin("cray-pals", cray_pals_ver))

-- Needed at runtime for nco utilities.
udunits_ver=os.getenv("udunits_ver") or "2.2.28"
load(pathJoin("udunits", udunits_ver))

-- Needed at runtime for nco utilities.
gsl_ver=os.getenv("gsl_ver") or "2.7"
load(pathJoin("gsl", gsl_ver))

nco_ver=os.getenv("nco_ver") or "5.0.6"
load(pathJoin("nco", nco_ver))

esmf_ver=os.getenv("esmf_ver") or "8.6.0"
load(pathJoin("esmf-C", esmf_ver))

whatis("Description: UFS_UTILS build environment")
