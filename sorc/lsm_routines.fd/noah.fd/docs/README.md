These directory contains copies of the noah LSM code needed to make any land state adjusments made neceessary by the DA updates. 
For noah, routines not labelled public/private, will shorty be replaced, and are very unlikely to change. For now, have made a second copy of the relevant routines. 

The copied versions are from:  https://github.com/NCAR/ccpp-physics -b master, commit 08b72bc1c23c48a823626d81f8e0a398685a35a3 (dated May, 2021). 
(used this version, as is formatted for doxygen. output zero diff to version in the global_workflow). 

Files directly copied:
namelist_soilveg.f  
set_soilveg.f   -> needed some changes for doxygen

sflx_snippet.f is a snippet of the file sflx.F, with minor changes to compil, and to enable frh2o routine to be called externally.

For complete documentation, see https://ufs-community.github.io/UFS_UTILS
