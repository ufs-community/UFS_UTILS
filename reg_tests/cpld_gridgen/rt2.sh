#!/bin/bash

set -x

readonly program=$(basename $0)
# PATHRT - Path to regression tests directory
readonly PATHRT="$(cd $(dirname $0) && pwd -P)"
export PATHRT
# PATHTR - Path to the UFS UTILS directory
readonly PATHTR="$(cd $PATHRT/../.. && pwd)"
export PATHTR

export compiler=${compiler:-intelllvm}
source $PATHTR/sorc/machine-setup.sh >/dev/null 2>&1
if [[ "$compiler" == "intelllvm" ]]; then
  if [[ ! -f ${PATHTR}/modulefiles/build.$target.$compiler.lua ]];then
     echo "IntelLLVM not available. Will use Intel Classic."
    compiler=intel
  fi
fi
echo "Machine: $target"
echo "Compiler: $compiler"

module use $PATHTR/modulefiles
module load build.$target.$compiler
if [[ $target = wcoss2 ]]; then
  module load netcdf
  module load nccmp
fi
set +x
module list
set -x

export CREATE_BASELINE=false
if [[ $target = hera ]]; then
  export MOM6_FIXDIR=/scratch1/NCEPDEV/global/glopara/fix/mom6/20220805
  STMP=${STMP:-/scratch2/NCEPDEV/stmp1/$USER}
  export NCCMP=nccmp
  BASELINE_ROOT=/scratch1/NCEPDEV/nems/role.ufsutils/ufs_utils/reg_tests/cpld_gridgen/baseline_data
fi


RUNDIR_ROOT=$STMP/CPLD_GRIDGEN/



declare -A tests
all_tests=""

i=0
while read -r line || [ "$line" ]; do

  line="${line#"${line%%[![:space:]]*}"}"
  [[ ${#line} == 0 ]] && continue
  [[ $line =~ \# ]] && continue

  TEST_NAME=$(echo $line | cut -d'|' -f1 | sed -e 's/^ *//' -e 's/ *$//')
  TEST_NAME=${TEST_NAME##mx}

  RUNDIR=$RUNDIR_ROOT/$TEST_NAME
  rm -fr $RUNDIR
  mkdir -p $RUNDIR
  export RUNDIR
  export OUTDIR_PATH=$RUNDIR
  export BASELINE=$BASELINE_ROOT/$TEST_NAME
  export REGRESSIONTEST_LOG=RegressionTests_$target.$compiler.${TEST_NAME}.log

  cp $PATHRT/parm/grid.nml.IN $RUNDIR
  cp $PATHTR/exec/cpld_gridgen $RUNDIR
  
  tests[$i]=$(sbatch --parsable --ntasks-per-node=1 --nodes=1 -t 0:10:00 -A fv3-cpu -q batch -J "test${i}" \
            -o log${i} -e log${i} $PATHTR/ush/cpld_gridgen.sh "$TEST_NAME")

  exit

  all_tests=${all_tests}":"${tests[$i]}

  ((i=i+1))
done < ./rt.conf

sbatch --nodes=1 -t 0:01:00 -A fv3-cpu -J summary -o logx -e logx \
       --open-mode=append -q batch -d afterok${all_tests} << EOF
#!/bin/bash
grep -a 'finished test' log*  > summary.log
EOF

exit
