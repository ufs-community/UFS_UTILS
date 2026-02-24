help([[ 
Load environment to compile UFS_UTILS on AWS EC2 using Intel
]])

prepend_path("MODULEPATH", "/opt/spack-stack/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")
prepend_path("MODULEPATH", "/opt/modulefiles")

stack_oneapi_ver=os.getenv("stack_oneapi_ver") or "2024.2.1"
stack_impi_ver=os.getenv("stack_impi_ver") or "2021.13"
cmake_ver=os.getenv("cmake_ver") or "3.27.9"

load(pathJoin("stack-oneapi", stack_oneapi_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))
load(pathJoin("cmake", cmake_ver))

local ufs_utils_modules = {
  {["jasper"]          = "2.0.32" },
  {["g2"]              = "3.5.1"  },
  {["ip"]              = "5.1.0" },
  {["sp"]              = "2.5.0" },
  {["netcdf-c"]        = "4.9.2"  },
  {["netcdf-fortran"]  = "4.6.1"  },
  {["bacio"]           = "2.4.1"  },
  {["nemsio"]          = "2.5.4"  },
  {["w3emc"]           = "2.10.0" },
  {["sigio"]           = "2.3.3"  },
  {["zlib-ng"]         = "2.2.1" },
  {["libpng"]          = "1.6.37" },
  {["hdf5"]            = "1.14.3" },
  {["nccmp"]           = "1.9.0.1"  },
  {["esmf"]            = "8.8.0"  },
  {["w3nco"]           = "2.4.1"  },
}

for i = 1, #ufs_utils_modules do
  for name, default_version in pairs(ufs_utils_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end

setenv("CC", "mpiicc")
setenv("CXX", "mpiicpc")
setenv("FC", "mpiifort")
