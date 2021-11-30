help([[
Load environment to compile UFS_UTILS emcsfc_ice_blend on Orion
]])

intel_ver=os.getenv("intel_ver") or "2020"
load(pathJoin("intel", intel_ver))

prepend_path("MODULEPATH", "/apps/contrib/NCEPLIBS/orion/modulefiles")

w3nco_ver=os.getenv("w3nco_ver") or "2.1.0"
load(pathJoin("w3nco", w3nco_ver))

bacio_ver=os.getenv("bacio_ver") or "2.2.0"
load(pathJoin("bacio", bacio_ver))

g2_ver=os.getenv("g2_ver") or "3.1.1"
load(pathJoin("g2", g2_ver))

jasper_ver=os.getenv("jasper_ver") or "1.900.1"
load(pathJoin("jasper", jasper_ver))
setenv("JASPER_LIB","${JASPER_LIBRARY_DIRS}/libjasper.a")

libpng_ver=os.getenv("libpng_ver") or "1.2.44"
load(pathJoin("png", libpng_ver))

zlib_ver=os.getenv("zlib_ver") or "1.2.6"
load(pathJoin("z", zlib_ver))

setenv("FCOMP","ifort")
setenv("FFLAGS","-O0 -i4")

whatis("Description: UFS_UTILS emcsfc_ice_blend build environment")
