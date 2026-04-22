#!/bin/bash

set -eux

#--------------------------------------------------------------------------
# This script runs once all tests have completed. It checks for
# failed tests and creates some final log files.
#--------------------------------------------------------------------------

cd $PATHRT

#--------------------------------------------------------------------------
# If there are any failed tests, combine the log files into one log
# called 'fail_test'.
#--------------------------------------------------------------------------

rm -f failed_tests
FAIL_FILES="fail_test_*"
for file in $FAIL_FILES; do
    if [[ -f "$file" ]]; then
        cat "$file" >> failed_tests
        rm -f $file
    fi
done

#--------------------------------------------------------------------------
# Combine the run logs for each test into one log.
#--------------------------------------------------------------------------

LOG_FILE=RegressionTests_${MACHINE_ID}.${compiler}.log
rm -f $LOG_FILE
for file in RegressionTests_${MACHINE_ID}.${compiler}.*.log
do
  if [[ -f "$file" ]]; then
    cat "$file" >> $LOG_FILE
    rm -f $file
  fi
done

#--------------------------------------------------------------------------
# Summarize the results in summary.log.
#--------------------------------------------------------------------------

rm -f summary.log
if [[ -e "failed_tests" ]]; then
  echo | tee -a $LOG_FILE
  cat failed_tests >>$LOG_FILE
  cat failed_tests >>summary.log
else
  echo | tee -a $LOG_FILE
  echo "REGRESSION TEST WAS SUCCESSFUL" | tee -a $LOG_FILE
  echo "All tests passed" >>summary.log
fi

exit 0
