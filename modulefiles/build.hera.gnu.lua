help([[
Load environment to compile UFS_UTILS on Hera using Gnu 13.3
]])

hpss_ver=os.getenv("hpss_ver") or ""
load(pathJoin("hpss", hpss_ver))

prepend_path("MODULEPATH", "/scratch2/NCEPDEV/stmp1/role.epic/installs/gnu/modulefiles")
prepend_path("MODULEPATH", "/scratch2/NCEPDEV/stmp1/role.epic/installs/openmpi/modulefiles")
prepend_path("MODULEPATH", "/scratch2/NCEPDEV/stmp1/role.epic/spack-stack/spack-stack-1.6.0_gnu13/envs/fms-2024.01/install/modulefiles/Core")

stack_gnu_ver=os.getenv("stack_gnu_ver") or "13.3.0"
load(pathJoin("stack-gcc", stack_gnu_ver))

stack_openmpi_ver=os.getenv("stack_openmpi_ver") or "4.1.6"
load(pathJoin("stack-openmpi", stack_openmpi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

bacio_ver=os.getenv("bacio_ver") or "2.4.1"
load(pathJoin("bacio", bacio_ver))

g2_ver=os.getenv("g2_ver") or "3.5.1"
load(pathJoin("g2", g2_ver))

ip_ver=os.getenv("ip_ver") or "4.3.0"
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

sigio_ver=os.getenv("sigio_ver") or "2.3.2"
load(pathJoin("sigio", sigio_ver))

nccmp_ver=os.getenv("nccmp_ver") or "1.9.0.1"
load(pathJoin("nccmp", nccmp_ver))

esmf_ver=os.getenv("esmf_ver") or "8.6.0"
load(pathJoin("esmf", esmf_ver))

nco_ver=os.getenv("nco_ver") or "5.0.6"
load(pathJoin("nco", nco_ver))


prepend_path("CPPFLAGS", " -I/apps/slurm_hera/23.11.3/include/slurm"," ")
prepend_path("LD_LIBRARY_PATH", "/apps/slurm_hera/23.11.3/lib")
setenv("LD_PRELOAD", "/scratch2/NCEPDEV/stmp1/role.epic/installs/gnu/13.3.0/lib64/libstdc++.so.6")

setenv("CC", "mpicc")
setenv("CXX", "mpic++")
setenv("FC", "mpif90")

whatis("Description: UFS_UTILS build environment")
