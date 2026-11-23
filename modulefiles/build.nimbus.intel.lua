help([[
Load environment to compile UFS_UTILS on NIMBUS
]])


source_sh("bash", "/opt/intel/oneapi/setvars.sh")
source_sh("bash", "/opt/intel/oneapi/compiler/2023.2.1/env/vars.sh")
source_sh("bash", "/opt/intel/oneapi/mpi/2021.18/env/vars.sh")

--append_path("MODULEPATH","/lfs/work/alexander_richert/stack/spack/var/spack/environments/nco-sci-intel-2021.10.0-v2/modules_flat/linux-rocky9-x86_64/Core")
append_path("MODULEPATH","/lfs/work/alexander_richert/stack/spack-stack/envs/nco-sci-intel-2021.10.0/modules_flat/Core")

cmake_ver=os.getenv("cmake_ver") or "3.31.11"
load(pathJoin("cmake", cmake_ver))

--PrgEnv_intel_ver=os.getenv("PrgEnv_intel_ver") or "8.3.3"
--load(pathJoin("PrgEnv-intel", PrgEnv_intel_ver))

--craype_ver=os.getenv("craype_ver") or "2.7.17"
--load(pathJoin("craype", craype_ver))

--intel_ver=os.getenv("intel_ver") or "2021.10.0-leeskur"
--load(pathJoin("intel-oneapi-compilers-classic", intel_ver))

--intel_mpi_ver=os.getenv("intel_mpi_ver") or "2021.18-qpoapwa"
--load(pathJoin("intel-oneapi-mpi", intel_mpi_ver))

-- Need the cray library path for C MPI libraries
--local cray_lib_path=os.getenv("CRAY_LD_LIBRARY_PATH") or ""
--prepend_path("LD_LIBRARY_PATH", cray_lib_path)

libjpeg_ver=os.getenv("libjpeg_ver") or "3.1.3"
load(pathJoin("libjpeg", libjpeg_ver))

zlib_ver=os.getenv("zlib_ver") or "1.3.2"
load(pathJoin("zlib", zlib_ver))

libpng_ver=os.getenv("libpng_ver") or "1.6.55"
load(pathJoin("libpng", libpng_ver))

--load("pkgconf/2.5.1-7zdb5kg")

hdf5_ver=os.getenv("hdf5_ver") or "1.14.5"
load(pathJoin("hdf5", hdf5_ver))

netcdf_c_ver=os.getenv("netcdf_c_ver") or "4.9.2"
load(pathJoin("netcdf-c", netcdf_c_ver))

netcdf_fortran_ver=os.getenv("netcdf_fortran_ver") or "4.6.1"
load(pathJoin("netcdf-fortran", netcdf_fortran_ver))

--netcdf_ver=os.getenv("pnetcdf_ver") or "1.12.2"
--load(pathJoin("pnetcdf-D", netcdf_ver))

bacio_ver=os.getenv("bacio_ver") or "2.4.1"
load(pathJoin("bacio", bacio_ver))

load("parallelio")
load("jasper")

-- Uncomment when CHGRES_ALL is ON.
--sfcio_ver=os.getenv("sfcio_ver") or "1.4.1"
--load(pathJoin("sfcio", sfcio_ver))

w3emc_ver=os.getenv("w3emc_ver") or "2.13.0"
load(pathJoin("w3emc", w3emc_ver))

nemsio_ver=os.getenv("nemsio_ver") or "2.5.5"
load(pathJoin("nemsio", nemsio_ver))

sigio_ver=os.getenv("sigio_ver") or "2.3.3"
load(pathJoin("sigio", sigio_ver))

sp_ver=os.getenv("sp_ver") or "2.5.0"
load(pathJoin("sp", sp_ver))

ip_ver=os.getenv("ip_ver") or "5.4.0"
load(pathJoin("ip", ip_ver))

--load("g2c/2.3.0-y5ymirz")

g2_ver=os.getenv("g2_ver") or "3.5.1"
load(pathJoin("g2", g2_ver))

-- Needed for mpiexec command.
--cray_pals_ver=os.getenv("cray_pals_ver") or "1.2.2"
--load(pathJoin("cray-pals", cray_pals_ver))

-- Needed at runtime for nco utilities.
udunits_ver=os.getenv("udunits_ver") or "2.2.28"
load(pathJoin("udunits", udunits_ver))

-- Needed at runtime for nco utilities.
gsl_ver=os.getenv("gsl_ver") or "2.8"
load(pathJoin("gsl", gsl_ver))

nco_ver=os.getenv("nco_ver") or "5.3.9"
load(pathJoin("nco", nco_ver))

load("python/3.11.15")
load("python-venv/1.0")
load("py-pyyaml/6.0.3")

esmf_ver=os.getenv("esmf_ver") or "8.8.0"
load(pathJoin("esmf", esmf_ver))

--nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
--load(pathJoin("nccmp-D", nccmp_ver))

load("openblas/0.3.33")

--setenv("CC", "icx")
--setenv("CXX", "icpx")
--setenv("FC", "ifx")
setenv("CC", "mpiicc")
setenv("CXX", "mpiicpx")
setenv("FC", "mpiifx")
setenv("CMAKE_Platform", "nimbus")

whatis("Description: UFS_UTILS build environment")
