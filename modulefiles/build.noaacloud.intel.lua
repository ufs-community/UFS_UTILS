help([[ 
Load environment to compile UFS_UTILS on NOAA CSPs using Intel
]])

prepend_path("MODULEPATH", "/contrib/spack-stack-rocky8/spack-stack-1.6.0/envs/ue-intel/install/modulefiles/Core")
prepend_path("MODULEPATH", "/apps/modules/modulefiles")
load("gnu")
load("stack-intel")
load("stack-intel-oneapi-mpi")
unload("gnu")
load("cmake/3.23.1")

load("common4noaacloud")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")

whatis("Description: UFS_UTILS build environment on NOAA Cloud")
