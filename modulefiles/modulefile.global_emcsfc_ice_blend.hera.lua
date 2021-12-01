help([[
Load environment to compile UFS_UTILS emcsfc_ice_blend on Hera
]])

intel_ver=os.getenv("intel_ver") or "18.0.5.274"
load(pathJoin("intel", intel_ver))

prepend_path("MODULEPATH", "/scratch2/NCEPDEV/nwprod/NCEPLIBS/modulefiles")

w3nco_ver=os.getenv("w3nco_ver") or "2.0.6"
load(pathJoin("w3nco", w3nco_ver))

bacio_ver=os.getenv("bacio_ver") or "2.0.2"
load(pathJoin("bacio", bacio_ver))

g2_ver=os.getenv("g2_ver") or "2.5.0"
load(pathJoin("g2", g2_ver))

jasper_ver=os.getenv("jasper_ver") or "1.900.1"
load(pathJoin("jasper", jasper_ver))

libpng_ver=os.getenv("libpng_ver") or "1.2.44"
load(pathJoin("png", libpng_ver))

zlib_ver=os.getenv("zlib_ver") or "1.2.11"
load(pathJoin("z", zlib_ver))

setenv("FCOMP","ifort")
setenv("FFLAGS","-O0 -i4")

whatis("Description: UFS_UTILS emcsfc_ice_blend build environment")
