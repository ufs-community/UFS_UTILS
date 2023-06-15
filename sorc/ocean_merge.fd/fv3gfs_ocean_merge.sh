#!/bin/bash


cat << EOF > input.nml

&mask_nml
 ocean_mask_dir="/scratch1/NCEPDEV/stmp4/Sanath.Kumar/CPLD_GRIDGEN/${ocn}/"
 ocnres="mx${ocn}"
 lake_mask_dir="/scratch2/NCEPDEV/stmp1/${USER}/my_grids/C${res}/"
 atmres="C${res}"
 out_dir= "/scratch2/NCEPDEV/stmp1/${USER}/ocean_merged/"
/  
EOF
echo "testing"
./merge

