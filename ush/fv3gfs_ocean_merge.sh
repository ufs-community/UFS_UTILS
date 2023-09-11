#!/bin/bash


    results_dir=$TEMP_DIR/ocean_merged/C${res}.mx${ocn}
    mkdir -p ${results_dir}

    cat << EOF > input.nml
     &mask_nml
     ocean_mask_dir="$ocean_mask_dir/${ocn}/"
     ocnres="mx${ocn}"
     lake_mask_dir="${TEMP_DIR}/C${res}/orog/"
     atmres="C${res}"
     out_dir="${results_dir}/"
    /
EOF


    time  ${exec_dir}/./ocean_merge
    set -x
    pwd

    echo run orog 2nd time

    for tnum in '1' '2' '3' '4' '5' '6'
    do
    cd ${TEMP_DIR}/C${res}/orog/tile$tnum
    echo $tnum $res $res 0 0 0 0 0 0 > INPS
    echo C${res}_grid.tile${tnum}.nc >> INPS

    echo none >> INPS
    echo ".false." >> INPS
    echo '"'${TEMP_DIR}/ocean_merged/C${res}.mx${ocn}/C${res}.mx${ocn}.tile${tnum}.nc'"' >> INPS

    cat INPS
    time ${exec_dir}/orog < INPS

    ncks -A -v lake_frac,lake_depth ${TEMP_DIR}/ocean_merged/C${res}.mx${ocn}/C${res}.mx${ocn}.tile${tnum}.nc out.oro.nc
    #cp out.oro.nc $out_dir/oro_C${res}.mx${ocn}.tile${tnum}.nc
    cp out.oro.nc $orog_dir/oro.C${res}.tile${tnum}.nc
    cp C${res}_grid.tile${tnum}.nc $out_dir/C${res}_grid.tile${tnum}.nc
    done
