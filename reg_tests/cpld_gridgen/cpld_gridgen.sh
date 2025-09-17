#!/bin/bash
set -eux

SECONDS=0

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
	-e "s/DO_POSTWGTS/$DO_POSTWGTS/g" \
	-e "s/ATMRESLIST/$ATMRESLIST/g"
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

        if [[ "$UPDATE_BASELINE" == "TRUE" ]] && [[ "${test_status}" == "FAIL" ]]; then
            source ${PATHTR}/reg_tests/get_hash.sh
            ${PATHTR}/reg_tests/update_baseline.sh ${BASELINE_ROOT}/.. "${TEST_NAME}" $commit_num
        fi

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

set +x
echo BEGIN cpld_gridgen.sh
set -x

cd $RUNDIR

RESNAME=${RESNAME:-$1}
TEST_NAME=$RESNAME
ATMLIST=${ATMLIST:-$2}
DEBUG=.false.
MASKEDIT=.false.
DO_POSTWGTS=.true.
MOSAICDIR_PATH=${MOSAICDIR_PATH:-$PATHTR/fix/orog}
FIXDIR_PATH=${MOM6_FIXDIR}/${RESNAME}
if [[ ${ATMLIST} -eq -1 ]]; then
   ATMRESLIST=12,24,48,96,192,384,768,1152,3072
else
   ATMRESLIST=${ATMLIST}
fi

APRUN=${APRUN:-"srun"}

if [ $RESNAME = 900 ]; then
    NI=40
    NJ=20
    TOPOGFILE=topog.nc
    EDITSFILE='none'
    if [ $DO_POSTWGTS == .true. ]; then
        #pre-generate SCRIP files for dst rectilinear grids using NCO
        $APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.9p00_SCRIP.nc -G latlon=20,40#lon_typ=grn_ctr#lat_typ=cap
    fi
fi

if [ $RESNAME = 500 ]; then
    NI=72
    NJ=35
    TOPOGFILE=ocean_topog.nc
    EDITSFILE='none'
    if [ $DO_POSTWGTS == .true. ]; then
	#pre-generate SCRIP files for dst rectilinear grids using NCO
        $APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.9p00_SCRIP.nc -G latlon=20,40#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.5p00_SCRIP.nc -G latlon=36,72#lon_typ=grn_ctr#lat_typ=cap
    fi
fi

if [ $RESNAME = 100 ]; then
    NI=360
    NJ=320
    MASKEDIT=.T.
    TOPOGFILE=topog.nc
    EDITSFILE=topo_edits_011818.nc
    if [ $DO_POSTWGTS == .true. ]; then
	#pre-generate SCRIP files for dst rectilinear grids using NCO
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.9p00_SCRIP.nc -G latlon=20,40#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.5p00_SCRIP.nc -G latlon=36,72#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.1p00_SCRIP.nc -G latlon=181,360#lon_typ=grn_ctr#lat_typ=cap
    fi
fi

if [ $RESNAME = 050 ]; then
    NI=720
    NJ=576
    TOPOGFILE=ocean_topog.nc
    EDITSFILE='none'
    if [ $DO_POSTWGTS == .true. ]; then
	#pre-generate SCRIP files for dst rectilinear grids using NCO
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.9p00_SCRIP.nc -G latlon=20,40#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.5p00_SCRIP.nc -G latlon=36,72#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.1p00_SCRIP.nc -G latlon=181,360#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.0p50_SCRIP.nc -G latlon=361,720#lon_typ=grn_ctr#lat_typ=cap
    fi
fi

if [ $RESNAME = 025 ]; then
    NI=1440
    NJ=1080
    TOPOGFILE=ocean_topog.nc
    EDITSFILE=All_edits.nc
    if [ $DO_POSTWGTS == .true. ]; then
	#pre-generate SCRIP files for dst rectilinear grids using NCO
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.9p00_SCRIP.nc -G latlon=20,40#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.5p00_SCRIP.nc -G latlon=36,72#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.1p00_SCRIP.nc -G latlon=181,360#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.0p50_SCRIP.nc -G latlon=361,720#lon_typ=grn_ctr#lat_typ=cap
	$APRUN -n 1 ncremap -g ${OUTDIR_PATH}/rect.0p25_SCRIP.nc -G latlon=721,1440#lon_typ=grn_ctr#lat_typ=cap
    fi
fi

edit_namelist < grid.nml.IN > grid.nml
$APRUN ./cpld_gridgen

# generate ice mesh
FSRC=${OUTDIR_PATH}/Ct.mx${RESNAME}_SCRIP_land.nc
FDST=${OUTDIR_PATH}/mesh.mx${RESNAME}.nc
$APRUN -n 1 ESMF_Scrip2Unstruct ${FSRC} ${FDST} 0

# generate kmt file for CICE
FSRC=${OUTDIR_PATH}/grid_cice_NEMS_mx${RESNAME}.nc
FDST=${OUTDIR_PATH}/kmtu_cice_NEMS_mx${RESNAME}.nc
ncks -O -v kmt ${FSRC} ${FDST}

check_results

elapsed_time=$( printf '%02dh:%02dm:%02ds\n' $((SECONDS%86400/3600)) $((SECONDS%3600/60)) $((SECONDS%60)) )
echo "Elapsed time: ${elapsed_time}. Have a nice day!" >> $PATHRT/${REGRESSIONTEST_LOG}
set +x
echo "Elapsed time: ${elapsed_time}. Have a nice day!"
