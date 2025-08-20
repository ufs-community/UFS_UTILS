help([[
Load environment to compile UFS_UTILS in a container using Intel
]])

prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.9.2/envs/unified-env/install/modulefiles/Core")

stack_oneapi_ver=os.getenv("stack_oneapi_ver") or "2024.2.0"
stack_impi_ver=os.getenv("stack_impi_ver") or "2021.13"

load(pathJoin("stack-oneapi", stack_oneapi_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

local ufs_utils_modules = {
  {["cmake"]           = "3.27.9" },
  {["bacio"]           = "2.4.1"  },
  {["g2"]              = "3.5.1"  },
  {["ip"]              = "5.1.0"  },
  {["sp"]              = "2.5.0"  },
  {["w3emc"]           = "2.10.0" },
  {["nemsio"]          = "2.5.4"  },
  {["sigio"]           = "2.3.3"  },
  {["libpng"]          = "1.6.37" },
  {["netcdf-c"]        = "4.9.2"  },
  {["netcdf-fortran"]  = "4.6.1"  },
  {["nccmp"]           = "1.9.0.1"},
  {["esmf"]            = "8.8.0"  },
  {["nco"]             = "5.2.4"  },
}

for i = 1, #ufs_utils_modules do
  for name, default_version in pairs(ufs_utils_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end

whatis("Description: UFS_UTILS build environment")
