#!/bin/bash
# A script to run jobs on flavor02

#nohup ./Projects/edmonds/IPASims/scripts/submit_flavor02_par-scans.sh > output.txt 2>&1 &

N_PARTICLES=2 # need to be one less so that seq works properly
PARTICLES=( "conversion" "proton" "deuteron" )
N_JOBS=( 5 5 5 );
N_EVTS=( 50000 6000000 6000000 );
#N_JOBS=( 2 2 2 );
#N_EVTS=( 100 100 100 );


TOTAL_JOBS=0
for I_PARTICLE in `seq 0 1 $N_PARTICLES`; do
    NAME="IPA-None_"${PARTICLES[$I_PARTICLE]}
    OUTDIR="data/batch/cd3-sims/"$NAME
    echo $NAME
    if [ -d $OUTDIR ] ; then
	echo "Already done "$NAME
	continue;
    else
	mkdir -p $OUTDIR
    fi
		
    for JOB_ID in `seq 1 1 ${N_JOBS[$I_PARTICLE]}`; do
	FILE_ID="Job"$JOB_ID"_"`date +%Y%m%d`"_"`date +%H%M%S`
		    
	logfile=$OUTDIR/"output_"$FILE_ID".log"
		    
	new_fclfile=$OUTDIR"/ipa-sims-"${PARTICLES[$I_PARTICLE]}"_"$FILE_ID".fcl"
	cat ~/Mu2e/Offline/Projects/edmonds-ipa-and-beam-test-sims/fcl/ipa-sims-${PARTICLES[$I_PARTICLE]}.fcl > $new_fclfile
	seed=$RANDOM
	echo -e "\nservices.user.SeedService.baseSeed         :  "$seed >> $new_fclfile
		    
	new_rootfile=$OUTDIR"/ipa-sims-"${PARTICLES[$I_PARTICLE]}"_"$FILE_ID".root"
	echo -e "\nservices.TFileService.fileName  : \""$new_rootfile"\"" >> $new_fclfile

	let TOTAL_JOBS=$TOTAL_JOBS+1;
	echo $new_rootfile

	nohup mu2e --config $new_fclfile --nevts ${N_EVTS[$I_PARTICLE]} > $logfile 2>&1 &

    done # jobs per particle

    while true ; do
	TEST=`ps -u aedmonds | grep mu2e | wc -l`
	if [ $TEST -eq 0 ] ; then
	    break;
	else
#		    echo "Waiting..."
	    sleep 10s
	fi
    done # while	    
done # particles

echo "Done!"

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
