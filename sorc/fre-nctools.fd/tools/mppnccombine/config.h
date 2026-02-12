#ifndef CONFIG_H
#define CONFIG_H

#define PACKAGE_NAME "FRE NCTools"
#define PACKAGE_VERSION "2024.05.02"
#define PACKAGE_BUGREPORT "https://github.com/NOAA-GFDL/FRE-NCtools/issues"

#define GIT_REVISION "unknown"
#define GIT_HEADHASH "unknown"
#define COPYRIGHT_YEAR "2024"

#define HAVE_GETRUSAGE 1

/* CRITICAL for 10GB files (Manual Addition) */
#define _FILE_OFFSET_BITS 64
#define _LARGEFILE64_SOURCE 1

#endif
