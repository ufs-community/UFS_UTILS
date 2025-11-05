help([[
Load environment to compile UFS_UTILS on Derecho using Intel Classic compilers
]])

setenv("LMOD_TMOD_FIND_FIRST","yes")
prepend_path("MODULEPATH", "/lustre/desc1/scratch/epicufsrt/contrib/modulefiles_extra")
prepend_path("MODULEPATH", "/glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_oneapi_ver") or "2024.2.1"
load(pathJoin("stack-oneapi", stack_oneapi_ver))

-- stack_oneapi_mpi_ver=os.getenv("stack_oneapi_mpi_ver") or "2021.13"
-- load(pathJoin("stack-intel-oneapi-mpi", stack_oneapi_mpi_ver))

stack_impi_ver=os.getenv("stack_cray_mpich_ver") or "8.1.29"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

bacio_ver=os.getenv("bacio_ver") or "2.4.1"
load(pathJoin("bacio", bacio_ver))

g2_ver=os.getenv("g2_ver") or "3.5.1"
load(pathJoin("g2", g2_ver))

ip_ver=os.getenv("ip_ver") or "5.1.0"
load(pathJoin("ip", ip_ver))

nemsio_ver=os.getenv("nemsio_ver") or "2.5.4"
load(pathJoin("nemsio", nemsio_ver))

sp_ver=os.getenv("sp_ver") or "2.5.0"
load(pathJoin("sp", sp_ver))

w3emc_ver=os.getenv("w3emc_ver") or "2.10.0"
load(pathJoin("w3emc", w3emc_ver))

-- Uncomment when CHGRES_ALL is ON
--sfcio_ver=os.getenv("sfcio_ver") or "1.4.1"
--load(pathJoin("sfcio", sfcio_ver))

sigio_ver=os.getenv("sigio_ver") or "2.3.3"
load(pathJoin("sigio", sigio_ver))

netcdf_c_ver=os.getenv("netcdf_c_ver") or "4.9.2"
load(pathJoin("netcdf-c", netcdf_c_ver))

netcdf_fortran_ver=os.getenv("netcdf_fortran_ver") or "4.6.1"
load(pathJoin("netcdf-fortran", netcdf_fortran_ver))

nccmp_ver=os.getenv("nccmp_ver") or "1.9.1.0"
load(pathJoin("nccmp", nccmp_ver))

esmf_ver=os.getenv("esmf_ver") or "8.6.1"
load(pathJoin("esmf", esmf_ver))

nco_ver=os.getenv("nco_ver") or "5.2.4"
load(pathJoin("nco", nco_ver))

-- setenv("I_MPI_CC", "icx")
-- setenv("I_MPI_FC", "ifort")

setenv("CC", "mpicc")
setenv("FC", "mpifort")

whatis("Description: UFS_UTILS build environment")
