help([[
Load environment to compile UFS_UTILS on Hera using Gnu
]])

hpss_ver=os.getenv("hpss_ver") or ""
load(pathJoin("hpss", hpss_ver))

prepend_path("MODULEPATH", "/scratch1/NCEPDEV/nems/role.epic/spack-stack/spack-stack-1.6.0/envs/unified-env/install/modulefiles/Core")
-- For openmpi:
prepend_path("MODULEPATH", "/scratch1/NCEPDEV/jcsda/jedipara/spack-stack/modulefiles")

stack_gcc_ver=os.getenv("stack_gcc_ver") or "9.2"
load(pathJoin("stack-gcc", gnu_ver))

stack_openmpi_ver=os.getenv("stack_openmpi_ver") or "4.1.5"
load(pathJoin("stack-openmpi", stack_openmpi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

bacio_ver=os.getenv("bacio_ver") or "2.4.1"
load(pathJoin("bacio", bacio_ver))

g2_ver=os.getenv("g2_ver") or "3.4.5"
load(pathJoin("g2", g2_ver))

ip_ver=os.getenv("ip_ver") or "4.3.0"
load(pathJoin("ip", ip_ver))

nemsio_ver=os.getenv("nemsio_ver") or "2.5.4"
load(pathJoin("nemsio", nemsio_ver))

sp_ver=os.getenv("sp_ver") or "2.5.0"
load(pathJoin("sp", sp_ver))

w3emc_ver=os.getenv("w3emc_ver") or "2.10.0"
load(pathJoin("w3emc", w3emc_ver))

sfcio_ver=os.getenv("sfcio_ver") or "1.4.1"
load(pathJoin("sfcio", sfcio_ver))

sigio_ver=os.getenv("sigio_ver") or "2.3.2"
load(pathJoin("sigio", sigio_ver))

hdf5_ver=os.getenv("hdf5_ver") or "1.14.0"
load(pathJoin("hdf5", hdf5_ver))

netcdf_c_ver=os.getenv("netcdf_c_ver") or "4.9.2"
load(pathJoin("netcdf-c", netcdf_c_ver))

netcdf_fortran_ver=os.getenv("netcdf_fortran_ver") or "4.6.1"
load(pathJoin("netcdf-fortran", netcdf_fortran_ver))

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

esmf_ver=os.getenv("esmf_ver") or "8.5.0"
load(pathJoin("esmf", esmf_ver))

nco_ver=os.getenv("nco_ver") or "5.0.6"
load(pathJoin("nco", nco_ver))

whatis("Description: UFS_UTILS build environment")
