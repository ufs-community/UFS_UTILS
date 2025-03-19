#!/bin/bash

set -eux

SECONDS=0

function edit_namelist {

    sed -e "s/SRCDIMS/$SRCDIMS/g" \
        -e "s/DSTDIMS/$DSTDIMS/g" \
	-e "s/FILETYPE/$FILETYPE/g" \
        -e "s|WTSDIR|$WTSDIR|g" \
	-e "s|GRDDIR|$GRDDIR|g" \
        -e "s/DO_DEBUG/$DO_DEBUG/g"
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

        for file in *mx*.nc; do
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

APRUN=${APRUN:-"srun --nodes=1 -A nems "}

# Two possible input files: ocean.nc and ice.nc
# One possible source grid, 1/4 deg (mx025: 1440,1080)
# Two possible destination grids, 1/2deg and 1deg (mx050: 720,576 and mx100: 360,320)
# The files produced will be ocean.mx[dest res].nc and ice.mx[des res].nc
# The program needs to execute twice, once for ocean and once for ice

# For the purposes of the RT, we can have staged input files retrieved from
# https://noaa-ufs-gefsv13replay-pds.s3.amazonaws.com/2021/03/2021032206/
# The required single mom file was created using
# ncks -O -v Temp,Salt,h,u replay-2021032206/MOM.res.nc ocean.nc
# ncks -v v,sfc -A replay-2021032206/MOM.res_1.nc ocean.nc
# I am assuming that the g-w will do the filename globbing and retrieval,
# process the NCO command and rename (timestamp) the output files at the end.

SRCDIMS="1440,1080"
FILETYPE=$FTYPE
WTSDIR=$WEIGHTS
GRDDIR=$WEIGHTS
DO_DEBUG=".false."

if [ $RESNAME = 050 ]; then
    DSTDIMS="720,576"
fi
if [ $RESNAME = 100 ]; then
    DSTDIMS="360,320"
fi

edit_namelist < ocniceprep.nml.IN > ocniceprep.nml

$APRUN ./oiprep

check_results

elapsed_time=$( printf '%02dh:%02dm:%02ds\n' $((SECONDS%86400/3600)) $((SECONDS%3600/60)) $((SECONDS%60)) )
echo "Elapsed time: ${elapsed_time}. Have a nice day!" >> $PATHRT/${REGRESSIONTEST_LOG}
set +x
echo "Elapsed time: ${elapsed_time}. Have a nice day!"

exit
