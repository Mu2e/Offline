#!/bin/bash
# A script to run jobs on flavor02

nohup ./Projects/edmonds/IPASims/scripts/submit_flavor02_par-scans.sh > output.txt 2>&1 &

#NJOBS=5;
#NEVTS=6000000;
#PARTICLE=Deuterons
#NAME=IPA-TDR_$PARTICLE
#OUTDIR="data/batch/"$NAME

#for job_id in `seq 1 1 $NJOBS`; do
#
 #   if [ ! -d $OUTDIR ] ; then
#	mkdir -p $OUTDIR
 #   fi
#
 #   file_id="Job"$job_id"_"`date +%Y%m%d`"_"`date +%H%M%S`
  #  echo "Starting run with id "$file_id
#
 #   logfile=$OUTDIR/"output_"$file_id".log"
#
 #   new_fclfile=$OUTDIR"/ProtonAbsorberRuns_"$PARTICLE"_"$file_id".fcl"
  #  cat ~/Mu2e/Offline/Projects/edmonds/IPASims/fcl/ProtonAbsorberRuns_$PARTICLE.fcl > $new_fclfile
   # seed=$RANDOM
    #echo -e "\nservices.user.SeedService.baseSeed         :  "$seed >> $new_fclfile
#
 #   new_rootfile=$OUTDIR"/PAbsRun_"$PARTICLE"_"$file_id".root"
  #  echo -e "\nservices.TFileService.fileName  : \""$new_rootfile"\"" >> $new_fclfile
#
#
 #   nohup mu2e --config $new_fclfile --nevts $NEVTS > $logfile 2>&1 &
#done
