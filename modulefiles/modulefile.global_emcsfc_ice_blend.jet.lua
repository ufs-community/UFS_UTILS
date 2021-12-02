help([[
Load environment to compile UFS_UTILS emcsfc_ice_blend on Jet
]])

intel_ver=os.getenv("intel_ver") or "18.0.5.274"
load(pathJoin("intel", intel_ver))

w3nco_ver=os.getenv("w3nco_ver") or "v2.0.6"
load(pathJoin("w3nco", w3nco_ver))

bacio_ver=os.getenv("bacio_ver") or "v2.0.2"
load(pathJoin("bacio", bacio_ver))

jasper_ver=os.getenv("jasper_ver") or "v1.900.1"
load(pathJoin("jasper", jasper_ver))

libpng_ver=os.getenv("libpng_ver") or "v1.2.44"
load(pathJoin("png", libpng_ver))

zlib_ver=os.getenv("zlib_ver") or "v1.2.6"
load(pathJoin("z", zlib_ver))

g2_ver=os.getenv("g2_ver") or "v3.1.0"
load(pathJoin("g2", g2_ver))

setenv("FCOMP","ifort")
setenv("FFLAGS","-O0 -i4")

whatis("Description: UFS_UTILS emcsfc_ice_blend build environment")
