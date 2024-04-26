#!/bin/bash
# 
# To generate the ocean mask 
#
# Check-out ufs-utils (and be sure to check out the ccpp submodule)
# cd fix 
# ./link_fixdirs.sh emc hera
# cd ../reg_tests/cpld_gridgen
# Edit the rt.conf and add the non-standard case(s) you want to generate
# Note you need to keep/run the C384_025 case because the lower resolution cases depend on it to generate the remapping weights used to create the CICE ICs and/or do the Post jobs.
# Edit rt.sh for proper accounts/partitions
# Build and run the test using ./rt.sh -b >output 2>&1 &
# Your results should be in /scratch1/NCEPDEV/stmp4/First.Last/CPLD_GRIDGEN/rt_#
#
#

    results_dir=$TEMP_DIR/ocean_merged/C${res}.mx${ocn}
    mkdir -p ${results_dir}
    
			
    cat << EOF > input.nml
     &mask_nml
     ocean_mask_dir="${home_dir}/fix/orog/C${res}/ocean_mask/${ocn}/"
     ocnres="mx${ocn}"
     lake_mask_dir="${TEMP_DIR}/C${res}/orog/"
     atmres="C${res}"
     out_dir="${results_dir}/"
     binary_lake=$binary_lake
    /
EOF

    echo Run ocean_merge program.
    time  ${exec_dir}/./ocean_merge
  
    rc=$?   

	if [[ $rc -ne 0 ]] ; then
     		echo FATAL ERROR running ocean_merge.
     		exit $rc
    	fi

    echo run orog 2nd time
    set +x
	
    for tnum in '1' '2' '3' '4' '5' '6'
    do
    cd ${TEMP_DIR}/C${res}/orog/tile$tnum

    echo C${res}_grid.tile${tnum}.nc > INPS
    echo ".false." >> INPS
    echo '"'${TEMP_DIR}/ocean_merged/C${res}.mx${ocn}/C${res}.mx${ocn}.tile${tnum}.nc'"' >> INPS

    cat INPS

    time ${exec_dir}/orog < INPS
    rc=$?   

    if [[ $rc -ne 0 ]] ; then
      echo "FATAL ERROR running orog."
      exit $rc
    fi

   ncks -4 -O ${TEMP_DIR}/ocean_merged/C${res}.mx${ocn}/C${res}.mx${ocn}.tile${tnum}.nc  ${TEMP_DIR}/ocean_merged/C${res}.mx${ocn}/C${res}.mx${ocn}.tile${tnum}.nc		
    ncks -A -v lake_frac,lake_depth ${TEMP_DIR}/ocean_merged/C${res}.mx${ocn}/C${res}.mx${ocn}.tile${tnum}.nc out.oro.nc
    cp out.oro.nc $orog_dir/oro.C${res}.tile${tnum}.nc
    done
