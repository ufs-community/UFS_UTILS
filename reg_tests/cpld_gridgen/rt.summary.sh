#!/bin/bash

set -x

cd $PATHRT

rm -f fail_test
FAIL_FILES="fail_test_*"
for file in $FAIL_FILES; do
    if [[ -f "$file" ]]; then
        cat "$file" >> fail_test
    fi
done

rm -f RegressionTests_$target.$compiler.log
for file in RegressionTests_$target.${compiler}.*.log
do
  if [[ -f "$file" ]]; then
    cat "$file" >> RegressionTests_$target.$compiler.log
    rm -f $file
  fi
done

rm -f summary.log
if [[ -e fail_test ]]; then
  echo | tee -a RegressionTests_$target.$compiler.log
    for file in fail_test_*; do
      cat $file >>RegressionTests_$target.$compiler.log
      cat $file >>summary.log
    done
else
  echo | tee -a RegressionTests_$target.$compiler.log
  echo "REGRESSION TEST WAS SUCCESSFUL" | tee -a RegressionTests_$target.$compiler.log
  echo "All tests passed" >>summary.log
fi
