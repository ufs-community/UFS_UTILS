help([[ 
Load environment to compile UFS_UTILS on NOAA CSPs using Intel
]])

prepend_path("MODULEPATH", "/contrib/spack-stack-rocky8/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")
prepend_path("MODULEPATH", "/apps/modules/modulefiles")

stack_oneapi_ver=os.getenv("stack_oneapi_ver") or "2024.2.1"
stack_impi_ver=os.getenv("stack_impi_ver") or "2021.13"
cmake_ver=os.getenv("cmake_ver") or "3.27.9"

load(pathJoin("stack-oneapi", stack_oneapi_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))
load(pathJoin("cmake", cmake_ver))

load("common4noaacloud")

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")

whatis("Description: UFS_UTILS build environment on NOAA Cloud")
