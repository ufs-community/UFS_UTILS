#!/bin/bash
set -eux

function edit_namelist {

    sed -e "s/NI_GLB/$NI/g" \
	-e "s/NJ_GLB/$NJ/g" \
	-e "s|FIXDIR|$FIXDIR_PATH|g" \
	-e "s|OUTDIR|$OUTDIR_PATH|g" \
	-e "s|MOSAICDIR|$MOSAICDIR_PATH|g" \
	-e "s/TOPOGFILE/$TOPOGFILE/g" \
	-e "s/EDITSFILE/$EDITSFILE/g" \
	-e "s/RESNAME/$RESNAME/g" \
	-e "s/DO_MASKEDIT/$MASKEDIT/g" \
	-e "s/DO_DEBUG/$DEBUG/g" \
	-e "s/DO_POSTWGTS/$DO_POSTWGTS/g"
}

check_results() {

    [ -o xtrace ] && set_x='set -x' || set_x='set +x'
    set +x

    local test_status=PASS
    # verification run
    if [[ $CREATE_BASELINE = false ]]; then

        echo | tee -a $PATHRT/$REGRESSIONTEST_LOG
        echo "Working dir = $RUNDIR" | tee -a $PATHRT/$REGRESSIONTEST_LOG
        echo "Baseline dir = $BASELINE" | tee -a $PATHRT/$REGRESSIONTEST_LOG
        echo | tee -a $PATHRT/$REGRESSIONTEST_LOG
        echo "Checking test $TEST_NAME results ...." | tee -a $PATHRT/$REGRESSIONTEST_LOG

        for file in $BASELINE/*.nc; do
            printf %s "Comparing " $(basename ${file}) "...." | tee -a $PATHRT/$REGRESSIONTEST_LOG

            if [[ ! -f $RUNDIR/$(basename ${file}) ]]; then
                echo "....MISSING file" | tee -a $PATHRT/$REGRESSIONTEST_LOG
                test_status=FAIL
            else
                $NCCMP -dmfqS -w format $(basename ${file}) $file >>${PATHRT}/nccmp_${TEST_NAME}.log 2>&1 && d=$? || d=$?
                if [[ $d -ne 0 ]]; then
                    echo "....NOT OK" | tee -a $PATHRT/$REGRESSIONTEST_LOG
                    test_status=FAIL
                else
                    echo "....OK" | tee -a $PATHRT/$REGRESSIONTEST_LOG
                fi
            fi
        done
        echo | tee -a $PATHRT/$REGRESSIONTEST_LOG
        # baseline creation run
    else

        echo | tee -a $PATHRT/$REGRESSIONTEST_LOG
        echo "Working dir = $RUNDIR" | tee -a $PATHRT/$REGRESSIONTEST_LOG
        echo "Moving baseline files to $NEW_BASELINE ...." | tee -a $PATHRT/$REGRESSIONTEST_LOG
        echo | tee -a $PATHRT/$REGRESSIONTEST_LOG

        mkdir -p $NEW_BASELINE

        for file in *.nc; do
            printf %s "Moving " $file "...." | tee -a $PATHRT/$REGRESSIONTEST_LOG

            cp $file $NEW_BASELINE/$file && d=$? || d=$?
            if [[ $d -ne 0 ]]; then
                echo "....NOT OK" | tee -a $PATHRT/$REGRESSIONTEST_LOG
                test_status=FAIL
            else
                echo "....OK" | tee -a $PATHRT/$REGRESSIONTEST_LOG
            fi
        done
        echo | tee -a $PATHRT/$REGRESSIONTEST_LOG

    fi

    if [[ $test_status == FAIL ]]; then
        echo "$TEST_NAME failed" >> $PATHRT/fail_test_$TEST_NAME
    fi
}

echo top of cpld_gridgen.sh

cd $RUNDIR

export RESNAME=${RESNAME:-$1}
TEST_NAME=$RESNAME
export DEBUG=.false.
export MASKEDIT=.false.
export DO_POSTWGTS=.true.
export MOSAICDIR_PATH=${MOSAICDIR_PATH:-$PATHTR/fix/orog}
export FIXDIR_PATH=${MOM6_FIXDIR}/${RESNAME}

APRUN=${APRUN:-"srun"}

if [ $RESNAME = 500 ]; then
    export NI=72
    export NJ=35
    export TOPOGFILE=ocean_topog.nc
    export EDITSFILE='none'
    if [ $DO_POSTWGTS == .true. ]; then
	#pre-generate SCRIP files for dst rectilinear grids using NCO
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.5p00_SCRIP.nc -G latlon=36,72#lon_typ=grn_ctr#lat_typ=cap
    fi
fi

if [ $RESNAME = 100 ]; then
    export NI=360
    export NJ=320
    export MASKEDIT=.T.
    export TOPOGFILE=topog.nc
    export EDITSFILE=topo_edits_011818.nc
    if [ $DO_POSTWGTS == .true. ]; then
	#pre-generate SCRIP files for dst rectilinear grids using NCO
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.5p00_SCRIP.nc -G latlon=36,72#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.1p00_SCRIP.nc -G latlon=181,360#lon_typ=grn_ctr#lat_typ=cap
    fi
fi

if [ $RESNAME = 050 ]; then
    export NI=720
    export NJ=576
    export TOPOGFILE=ocean_topog.nc
    export EDITSFILE='none'
    if [ $DO_POSTWGTS == .true. ]; then
	#pre-generate SCRIP files for dst rectilinear grids using NCO
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.5p00_SCRIP.nc -G latlon=36,72#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.1p00_SCRIP.nc -G latlon=181,360#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.0p50_SCRIP.nc -G latlon=361,720#lon_typ=grn_ctr#lat_typ=cap
    fi
fi

if [ $RESNAME = 025 ]; then
    export NI=1440
    export NJ=1080
    export TOPOGFILE=ocean_topog.nc
    export EDITSFILE=All_edits.nc
    if [ $DO_POSTWGTS == .true. ]; then
	#pre-generate SCRIP files for dst rectilinear grids using NCO
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.5p00_SCRIP.nc -G latlon=36,72#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.1p00_SCRIP.nc -G latlon=181,360#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.0p50_SCRIP.nc -G latlon=361,720#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.0p25_SCRIP.nc -G latlon=721,1440#lon_typ=grn_ctr#lat_typ=cap
    fi
fi

edit_namelist < grid.nml.IN > grid.nml
$APRUN ./cpld_gridgen

# generate ice mesh
export FSRC=${OUTDIR_PATH}/Ct.mx${RESNAME}_SCRIP_land.nc
export FDST=${OUTDIR_PATH}/mesh.mx${RESNAME}.nc
$APRUN -n 1 ESMF_Scrip2Unstruct ${FSRC} ${FDST} 0

# generate kmt file for CICE
export FSRC=${OUTDIR_PATH}/grid_cice_NEMS_mx${RESNAME}.nc
export FDST=${OUTDIR_PATH}/kmtu_cice_NEMS_mx${RESNAME}.nc
ncks -O -v kmt ${FSRC} ${FDST}

check_results
