whatis("Description: UFS_UTILS build environment common libraries")

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
  {["zlib"]            = "1.2.13" },
  {["libpng"]          = "1.6.37" },
  {["hdf5"]            = "1.14.3" },
  {["netcdf"]          = "4.7.0"  },
  {["nccmp"]           = "1.9.1"  },
  {["esmf"]            = "8.6.1"  },
  {["nco"]             = "5.1.6"  },
}

for i = 1, #ufs_utils_modules do
  for name, default_version in pairs(ufs_utils_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end
