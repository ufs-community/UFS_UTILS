#!/bin/sh

#------------------------------------------------------------------------
#
# Script name: fv3gfs_chgres.sh
#
# Description: Creates initial conditions for coldstarting the fv3
#              forecast model.
#
# How to use this script:
#
#  - Uncomment the job cards for your machine.  Ensure the
#    'machine' variable is uncommented.  Ensure job cards for
#    the other machines are prefixed by "##".  
#  - Check the job card account.  On Hera, this is 
#    "#SBATCH -q $account".  On WCOSS, this is "#BSUB -P $account".
#    You may not have permissions to run under the default account.
#  - Make sure you have linked the 'fixed' directories using
#    sorc/link_fixdirs.sh.
#  - Set the following shell variables:
#      PSLOT     - Your experiment name
#      CASE_HIGH - Deterministic resolution
#      CASE_ENKF - Enkf resolution
#      PTMP      - Directory were the initial condition files will be 
#                  placed.
#      CDUMP     - Set to 'gdas' for cycling with enkf.  Set to 'gfs'
#                  for free forecast only.  NOTE: only the 'gfs'
#                  cycle is supported when using fv3gfs data as
#                  input.
#      CDATE     - year/month/day/hour of your experiment.
#      HOMEgfs   - Path to your checkout repository.
#  - Start the script as follows:
#      Hera: "sbatch fv3gfs_chgres.sh"
#      WCOSS Dell and Cray: "cat fv3gfs_chgres.sh | bsub"
#
#------------------------------------------------------------------------

#----WCOSS_DELL JOBCARD
#BSUB -P FV3GFS-T2O
#BSUB -o log.chgres.%J
#BSUB -e log.chgres.%J
#BSUB -J fv3_chgres
#BSUB -q dev
#BSUB -M 6000
#BSUB -W 10:00
#BSUB -R span[ptile=14]
#BSUB -n 14
export machine=WCOSS_DELL_P3

#----WCOSS_CRAY JOBCARD
##BSUB -L /bin/sh
##BSUB -P FV3GFS-T2O
##BSUB -oo log.chgres.%J
##BSUB -eo log.chgres.%J
##BSUB -J fv3_chgres
##BSUB -q dev
##BSUB -M 2400
##BSUB -W 10:00
##BSUB -extsched 'CRAYLINUX[]'
#export machine=WCOSS_C

#---- Hera JOBCARD
#---- Submit as: sbatch $script
##SBATCH -J fv3_chgres_driver
##SBATCH -A fv3-cpu
##SBATCH --open-mode=truncate
##SBATCH -o log.chgres
##SBATCH -e log.chgres
##SBATCH --nodes=1
##SBATCH -q batch
##SBATCH -t 03:00:00
#export machine=HERA

set -x

export PSLOT=fv3test
export CDUMP=gfs
export CASE_HIGH=C768
export CASE_ENKF=C384
export CDATE=2019061000

if [ $machine = WCOSS_C ]; then
 export HOMEgfs=$LS_SUBCWD/..
 export PTMP=${PTMP:-"/gpfs/hps2/ptmp/$USER"}
 export zero_bias_dir=/gpfs/hps3/emc/global/noscrub/emc.glopara/ICS/bias_zero
elif [ $machine = WCOSS_DELL_P3 ]; then
 export HOMEgfs=$LS_SUBCWD/..
 export PTMP=${PTMP:-"/gpfs/dell2/ptmp/$USER"}
 export zero_bias_dir=/gpfs/hps3/emc/global/noscrub/emc.glopara/ICS/bias_zero
elif [ $machine = HERA ]; then
 export HOMEgfs=$SLURM_SUBMIT_DIR/..
 export PTMP=${PTMP:-"/scratch2/NCEPDEV/stmp1/$USER"}
 export zero_bias_dir=/scratch1/NCEPDEV/global/glopara/bias_zero
else 
 echo "Please define HOMEgfs and PTMP for your machine. exit"
 exit
fi

#-------------------------------------------------------------------------------------------------
# Perform simple check to ensure user has HOMEgfs set correctly and has linked the 'fix'
# directories.
#-------------------------------------------------------------------------------------------------

if [ ! -d ${HOMEgfs}/fix/fix_fv3/$CASE_HIGH ]; then
  set +x
  echo FATAL ERROR: CANT FIND FIX DIRECTORIES. 
  echo ENSURE HOMEgfs PATH IS CORRECT.
  echo YOUR HOMEgfs SETTING IS: $HOMEgfs
  echo ALSO ENSURE sorc/link_fixdirs.sh WAS RUN.
  exit 1
fi

export NSTSMTH=NO                                  ##apply 9-point smoothing to nsst tref
export NST_TF_CHG=$HOMEgfs/exec/nst_tf_chg.x
export ZERO_BIAS=YES                                ##zeroed out all bias and radsat files 


export ymd=`echo $CDATE | cut -c 1-8`
export cyc=`echo $CDATE | cut -c 9-10`
export yy=`echo $CDATE | cut -c 1-4`
export mm=`echo $CDATE | cut -c 5-6`
export dd=`echo $CDATE | cut -c 7-8`

export ROTDIR=$PTMP/$PSLOT
export RUNDIR=$ROTDIR/chgres

export NODES=1
export APRUNC=""

. $HOMEgfs/ush/load_fv3gfs_modules.sh
status=$?
[[ $status -ne 0 ]] && exit $status

if [ $machine = WCOSS_C ]; then
 export KMP_AFFINITY=disabled
 export OMP_NUM_THREADS_CH=24
 export APRUNC="aprun -n 1 -N 1 -j 1 -d $OMP_NUM_THREADS_CH -cc depth"
 export APRUNTF='aprun -q -j1 -n1 -N1 -d1 -cc depth'
 export SUB=$LS_SUBCWD/../util/sub_wcoss_c
 export ACCOUNT=FV3GFS-T2O
 export QUEUE=dev
 export QUEUE_TRANS=dev_transfer 
elif [ $machine = WCOSS_DELL_P3 ]; then
 export OMP_NUM_THREADS_CH=14
 export KMP_AFFINITY=disabled
 export APRUNC="time"
 export APRUNTF="time"
 export SUB=$LS_SUBCWD/../util/sub_wcoss_d
 export ACCOUNT=FV3GFS-T2O
 export QUEUE=dev
 export QUEUE_TRANS=dev_transfer 
elif [ $machine = HERA ]; then
 export APRUNC=time
 export APRUNTF=time
 export OMP_NUM_THREADS_CH=24
 export SUB=$SLURM_SUBMIT_DIR/../util/sub_slurm
 export ACCOUNT=fv3-cpu
 export QUEUE=batch
 export QUEUE_TRANS=batch
else
 echo "$machine not supported, exit"
 exit
fi

#----------------------------
#----------------------------
#--for high-res
#----------------------------
#----------------------------

export CASE=$CASE_HIGH
export INIDIR=$RUNDIR/$CDUMP/$CASE
export COMROT=$ROTDIR/${CDUMP}.$ymd/$cyc
export OUTDIR=$COMROT/INPUT
export DATA=$INIDIR/stmp
rm -rf $INIDIR $OUTDIR $DATA
mkdir -p $INIDIR $OUTDIR  $DATA
cd $INIDIR ||exit 8

#--------------------------------------------------------------------------
# See if requested data is in the ops directory.  If not, pull from hpss.
#--------------------------------------------------------------------------

if [ -s $COMROOT/gfs/prod/${CDUMP}.${ymd}/$cyc ]; then

   if [ $CDUMP = gdas ]; then
     set +x
     echo GDAS CYCLE NOT SUPPORTED.
     echo USE NEW PARALLEL VERSION OF CHGRES.
     echo CONTACT REPOSITORY MANAGER FOR DETAILS.
     exit
   fi
   atm=./${CDUMP}.t${cyc}z.atmanl.nemsio
   sfc=./${CDUMP}.t${cyc}z.sfcanl.nemsio
   nst=$sfc  # for fv3gfs, nst fields stored in surface file
   for ff in $atm $sfc ; do
     cp $COMROOT/gfs/prod/${CDUMP}.${ymd}/$cyc/$ff .
   done

else   ##get data from HPSS archive

  if [ $CDATE -lt 2014050100 ]; then

    HPSSPATH=/NCEPPROD/hpssprod/runhistory/rh$yy/$yy$mm/$yy$mm$dd      ##use operational nems gfs nems gfs ics
    nst=" "
    if [ $CDUMP = gdas ]; then
      tarball_high=com_gfs_prod_${CDUMP}.${CDATE}.tar
      atm=./gdas1.t${cyc}z.sanl
      sfc=./gdas1.t${cyc}z.sfcanl
      biascr=./gdas1.t${cyc}z.abias
      biascr_pc=" "
      radstat=./gdas1.t${cyc}z.radstat
      aircraft_t_bias=" "
    else
      tarball_high=com_gfs_prod_${CDUMP}.${CDATE}.anl.tar
      atm=./gfs.t${cyc}z.sanl
      sfc=./gfs.t${cyc}z.sfcanl
    fi

  elif [ $CDATE -le 2017072012 ]; then
    if [ $CDATE -ge 2016110100 ]; then
      oldexp=prnemsrn
    elif [ $CDATE -ge 2016050100 ]; then
      oldexp=pr4rn_1605
    elif [ $CDATE -ge 2015121500 ]; then
      oldexp=pr4rn_1512
    elif [ $CDATE -ge 2015050200 ]; then
      oldexp=pr4rn_1505
    elif [ $CDATE -ge 2014073000 ]; then
      oldexp=pr4rn_1408
    elif [ $CDATE -ge 2014050100 ]; then
      oldexp=pr4rn_1405
    else
      echo "NEMS GSM retro ICs do not exit, exit"
      exit 1
    fi

    HPSSPATH=/5year/NCEPDEV/emc-global/emc.glopara/WCOSS_C/$oldexp    ##use q3fy17 nems gfs parallel ics
    tarball_high=${CDATE}${CDUMP}.tar
    atm=gfnanl.${CDUMP}.$CDATE
    sfc=sfnanl.${CDUMP}.$CDATE
    nst=nsnanl.${CDUMP}.$CDATE
    biascr=biascr.${CDUMP}.$CDATE
    biascr_pc=biascr_pc.${CDUMP}.$CDATE
    aircraft_t_bias=aircraft_t_bias.${CDUMP}.$CDATE
    radstat=radstat.${CDUMP}.$CDATE
  elif [[ $CDATE -le 2019061118 ]]; then
    HPSSPATH=/NCEPPROD/hpssprod/runhistory/rh$yy/$yy$mm/$yy$mm$dd      ##use operational nems gfs nems gfs ics
    if [ $CDUMP = gfs ]; then
      tarball_high=gpfs_hps_nco_ops_com_gfs_prod_${CDUMP}.${CDATE}.anl.tar
    else
      tarball_high=gpfs_hps_nco_ops_com_gfs_prod_${CDUMP}.${CDATE}.tar
    fi
    atm=./${CDUMP}.t${cyc}z.atmanl.nemsio
    sfc=./${CDUMP}.t${cyc}z.sfcanl.nemsio
    nst=./${CDUMP}.t${cyc}z.nstanl.nemsio
    biascr=./${CDUMP}.t${cyc}z.abias 
    biascr_pc=./${CDUMP}.t${cyc}z.abias_pc
    aircraft_t_bias=./${CDUMP}.t${cyc}z.abias_air  
    radstat=./${CDUMP}.t${cyc}z.radstat
  else  # FV3GFS starting 2019061200
    HPSSPATH=/NCEPPROD/hpssprod/runhistory/rh$yy/$yy$mm/$yy$mm$dd
    if [ $CDUMP = gfs ]; then
      tarball_high=gpfs_dell1_nco_ops_com_gfs_prod_${CDUMP}.${yy}${mm}${dd}_${cyc}.${CDUMP}_nemsioa.tar
    else
      set +x
      echo GDAS CYCLE NOT SUPPORTED.
      echo USE NEW PARALLEL VERSION OF CHGRES.
      echo CONTACT REPOSITORY MANAGER FOR DETAILS.
      exit
    fi
    atm=./${CDUMP}.$yy$mm$dd/$cyc/${CDUMP}.t${cyc}z.atmanl.nemsio
    sfc=./${CDUMP}.$yy$mm$dd/$cyc/${CDUMP}.t${cyc}z.sfcanl.nemsio
  fi

#--extract ICs from hpss
cat > read_hpss.sh <<EOF1

   . $HOMEgfs/ush/load_fv3gfs_modules.sh   2>>/dev/null

   cd $INIDIR
   htar -xvf  $HPSSPATH/$tarball_high $atm $sfc $nst  
   if [[ $CDATE -ge 2019061200 ]]; then
     ln -fs $atm .
     ln -fs $sfc .
   fi
   if [ $CDUMP = gdas ]; then
      cd $COMROT 
      htar -xvf  $HPSSPATH/$tarball_high $biascr $biascr_pc $aircraft_t_bias $radstat

      [[ '$biascr' = ./gdas1.t${cyc}z.abias ]] && mv $biascr ./gdas.t${cyc}z.abias
      [[ '$radstat' = ./gdas1.t${cyc}z.radstat ]] && mv $radstat ./gdas.t${cyc}z.radstat

      [[ '$biascr' = biascr.${CDUMP}.$CDATE ]] && mv biascr.${CDUMP}.$CDATE  ${CDUMP}.t${cyc}z.abias 
      [[ '${biascr_pc}' = biascr_pc.${CDUMP}.$CDATE ]] &&  mv biascr_pc.${CDUMP}.$CDATE ${CDUMP}.t${cyc}z.abias_pc
      [[ '$aircraft_t_bias' = aircraft_t_bias.${CDUMP}.$CDATE ]] && mv aircraft_t_bias.${CDUMP}.$CDATE ${CDUMP}.t${cyc}z.abias_air
      [[ '$radstat' = radstat.${CDUMP}.$CDATE ]] &&  mv radstat.${CDUMP}.$CDATE  ${CDUMP}.t${cyc}z.radstat

   fi
EOF1
chmod u+x read_hpss.sh
$SUB -a $ACCOUNT -q $QUEUE_TRANS -p 1/1/S -r 1024/1/1 -t 2:00:00 -j read_hpss -o read_hpss.out $INIDIR/read_hpss.sh

#................................................
fi
#................................................

testfile=$INIDIR/$sfc                 
nsleep=0; tsleep=120;  msleep=50
while test ! -s $testfile -a $nsleep -lt $msleep;do
  sleep $tsleep; nsleep=`expr $nsleep + 1`
done
sleep 300

if [ ! -s $testfile ]; then 
  echo "$testfile does not exist, exit !"
  exit 1
fi

[[ $atm = ./gdas1.t${cyc}z.sanl ]] && mv $INIDIR/$atm $INIDIR/gdas.t${cyc}z.sanl
[[ $sfc = ./gdas1.t${cyc}z.sfcanl ]] && mv $INIDIR/$sfc $INIDIR/gdas.t${cyc}z.sfcanl

#------------------------------
if [ $NSTSMTH = "YES" ]; then
#------------------------------
mv $nst fnsti
rm -f tf_chg_parm.input
cat >tf_chg_parm.input <<EOF1
&config
nsmth=3,istyp=0,
/
EOF1
$APRUNTF $NST_TF_CHG <tf_chg_parm.input 
mv fnsto $nst
if [ $? -ne 0 ] ; then
  echo "NST_TF_CHG for $CDUMP $CASE failed. exit"
  exit 1
fi
#------------------------------
fi
#------------------------------

$HOMEgfs/ush/global_chgres_driver.sh
if [ $? -ne 0 ] ; then
 echo "chgres for $CDUMP $CASE failed. exit"
 exit 1
fi

if [ $CDUMP = "gdas" -a $ZERO_BIAS = "YES" ]; then
   cd $COMROT
   for ff in abias abias_air abias_pc radstat; do
    mv gdas.t${cyc}z.${ff}  gdas.t${cyc}z.${ff}_save
    cp -p ${zero_bias_dir}/gdas.t18z.${ff} gdas.t${cyc}z.${ff}
   done
fi 

[[ $CDUMP = gfs ]] && exit
[[ $CDATE -lt 2012052100 ]] && exit  # no enkf data before this date

#----------------------------
#----------------------------
#--for ENKF 
#----------------------------
#----------------------------
export CASE=$CASE_ENKF
export INIDIR=$RUNDIR/$CDUMP/$CASE
rm -rf $INIDIR; mkdir -p $INIDIR 
cd $INIDIR ||exit 8

if [ $CDATE -lt 2014050100 ]; then
  HPSSPATH=/NCEPPROD/hpssprod/runhistory/rh$yy/$yy$mm/$yy$mm$dd      ##use operational nems gfs nems gfs ics
  tarball_enkf_atm=com_gfs_prod_enkf.${ymd}_${cyc}.anl.tar    
  testfile=$INIDIR/siganl_${CDATE}_mem080
elif [ $CDATE -le 2017072012 ]; then
   if [ $CDATE -ge 2016110100 ]; then
     oldexp=prnemsrn
   elif [ $CDATE -ge 2016050100 ]; then
     oldexp=pr4rn_1605
   elif [ $CDATE -ge 2015121500 ]; then
     oldexp=pr4rn_1512
   elif [ $CDATE -ge 2015050200 ]; then
     oldexp=pr4rn_1505
   elif [ $CDATE -ge 2014073000 ]; then
     oldexp=pr4rn_1408
   elif [ $CDATE -ge 2014050100 ]; then
     oldexp=pr4rn_1405
   else
     echo "NEMS GSM retro ICs do not exit, exit"
     exit 1
   fi

 HPSSPATH=/5year/NCEPDEV/emc-global/emc.glopara/WCOSS_C/$oldexp    ##use q3fy17 nems gfs parallel ics
 tarball_enkf_atm=${CDATE}gdas.enkf.anl.tar
 tarball_enkf_sfcnst=${CDATE}gdas.enkf.sfcanl.tar
 testfile=$INIDIR/sfcanl_${CDATE}_mem080
else
 HPSSPATH=/NCEPPROD/hpssprod/runhistory/rh$yy/$yy$mm/$yy$mm$dd      ##use operational nems gfs nems gfs ics
 tarball_enkf_atm=gpfs_hps_nco_ops_com_gfs_prod_enkf.${ymd}_${cyc}.anl.tar    
 testfile=$INIDIR/gdas.t${cyc}z.sfcanl.mem080.nemsio
fi

#--extract ICs from hpss
cat > read_hpss.sh <<EOF

   . $HOMEgfs/ush/load_fv3gfs_modules.sh   2>>/dev/null

   cd $INIDIR
   htar -xvf  $HPSSPATH/$tarball_enkf_atm 
   if [ $CDATE -le 2017072000 ] && [ $CDATE -gt 2014050100 ]; then
      htar -xvf  $HPSSPATH/$tarball_enkf_sfcnst 
   fi
EOF
chmod u+x read_hpss.sh
$SUB -a $ACCOUNT -q $QUEUE_TRANS -p 1/1/S -r 1024/1/1 -t 2:00:00 -j read_hpss -o read_hpss.out $INIDIR/read_hpss.sh

nsleep=0; tsleep=120;  msleep=50
while test ! -s $testfile -a $nsleep -lt $msleep;do
  sleep $tsleep; nsleep=`expr $nsleep + 1`
done

if [ ! -s $testfile ]; then
 echo " $testfile does not exist, exit"
 exit
fi

#---------------------------
n=1
while [ $n -le 80 ]; do
#---------------------------
mem=$(printf %03i $n)
mchar=mem$(printf %03i $n)


export COMROT=$ROTDIR/enkfgdas.$ymd/$cyc/$mchar
export OUTDIR=$COMROT/INPUT
export DATA=$INIDIR/$mchar
rm -rf $OUTDIR $DATA
mkdir -p $OUTDIR $DATA

if [ $CDATE -lt 2014050100 ]; then
  atm=gdas.t${cyc}z.sanl
  sfc=gdas.t${cyc}z.sfcanl
  rm -fr $atm $sfc 
else
  atm=${CDUMP}.t${cyc}z.atmanl.nemsio
  sfc=${CDUMP}.t${cyc}z.sfcanl.nemsio
  nst=${CDUMP}.t${cyc}z.nstanl.nemsio
  rm -rf $atm $sfc $nst 
fi

if [ $CDATE -lt 2014050100 ]; then
   ln -fs siganl_${CDATE}_$mchar $atm
   ln -fs sfcanl_${CDATE}_$mchar $sfc
elif [ $CDATE -le 2017072000 ]; then
   ln -fs siganl_${CDATE}_$mchar $atm
   ln -fs sfcanl_${CDATE}_$mchar $sfc
   ln -fs nstanl_${CDATE}_$mchar $nst
   cp     nstanl_${CDATE}_$mchar fnsti 
else
   ln -fs gdas.t${cyc}z.ratmanl.${mchar}.nemsio $atm
   ln -fs gdas.t${cyc}z.sfcanl.${mchar}.nemsio $sfc
   ln -fs gdas.t${cyc}z.nstanl.${mchar}.nemsio $nst
   cp     gdas.t${cyc}z.nstanl.${mchar}.nemsio fnsti
fi

if [ $NSTSMTH = "YES" ]; then

rm -f tf_chg_parm.input
cat >tf_chg_parm.input <<EOF1
&config
nsmth=3,istyp=0,
/
EOF1
$APRUNTF $NST_TF_CHG <tf_chg_parm.input
rm -f $nst
mv fnsto $nst
if [ $? -ne 0 ] ; then
  echo "NST_TF_CHG for $CDUMP $CASE failed. exit"
  exit 1
fi

fi # NSTSMTH IS YES

$HOMEgfs/ush/global_chgres_driver.sh
if [ $? -ne 0 ] ; then
 echo "chgres for enkf $CASE ${mchar} failed. exit"
 exit
fi

#---------------------------
n=$((n+1))
done
#---------------------------

exit
