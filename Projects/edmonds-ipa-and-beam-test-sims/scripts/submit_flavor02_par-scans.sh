#!/bin/bash
# A script to run jobs on flavor02

RMIN=200;
RMAX=200;
RSTEP=100;

HALF_L_MIN=250;
HALF_L_MAX=500;
HALF_L_STEP=125;

ZMIN=500;
ZMAX=1000;
ZSTEP=250;

N_PARTICLES=2 # need to be one less so that seq works properly
PARTICLES=( "ConvEMinus" "Protons" "Deuterons" )
N_JOBS=( 1 2 2 );
N_EVTS=( 50000 3000000 3000000 );

ipa_geom_file="Mu2eG4/geom/protonAbsorber_cylindrical_v2.txt"

for R in `seq $RMIN $RSTEP $RMAX`; do
    for HALF_L in `seq $HALF_L_MIN $HALF_L_STEP $HALF_L_MAX`; do
	for Z in `seq $ZMIN $ZSTEP $ZMAX`; do
	    # Change the IPA geometry
	    # Z Position
	    val_to_change_to=$Z";"
	    line=`cat $ipa_geom_file | grep "protonabsorber.distFromTargetEnd"`
	    val_to_change=`echo $line | cut -d ' ' -f 4`
	    newline=${line/$val_to_change/$val_to_change_to} # change the parameter in the line
	    sed -i "s/$line/$newline/g" $ipa_geom_file # change the line in the file

	    # Radius (need to change both the front and back)
	    val_to_change_to=$R";"
	    line=`cat $ipa_geom_file | grep "protonabsorber.OutRadius0"`
	    val_to_change=`echo $line | cut -d ' ' -f 4`
	    newline=${line/$val_to_change/$val_to_change_to} # change the parameter in the line
	    sed -i "s/$line/$newline/g" $ipa_geom_file # change the line in the file

	    val_to_change_to=$R";"
	    line=`cat $ipa_geom_file | grep "protonabsorber.OutRadius1"`
	    val_to_change=`echo $line | cut -d ' ' -f 4`
	    newline=${line/$val_to_change/$val_to_change_to} # change the parameter in the line
	    sed -i "s/$line/$newline/g" $ipa_geom_file # change the line in the file

	    val_to_change_to=$HALF_L";"
	    line=`cat $ipa_geom_file | grep "protonabsorber.halfLength"`
	    val_to_change=`echo $line | cut -d ' ' -f 4`
	    newline=${line/$val_to_change/$val_to_change_to} # change the parameter in the line
	    sed -i "s/$line/$newline/g" $ipa_geom_file # change the line in the file


	    TOTAL_JOBS=0
	    for I_PARTICLE in `seq 0 1 $N_PARTICLES`; do
		NAME="IPA-Cyl_R"$R"mm_L"$HALF_L"mm_Z"$Z"mm_"${PARTICLES[$I_PARTICLE]}
		OUTDIR="data/batch/"$NAME
#		echo $NAME
		if [ -d $OUTDIR ] ; then
		    echo "Already done "$NAME
		    continue;
		else
		    mkdir -p $OUTDIR
		fi
		
		for JOB_ID in `seq 1 1 ${N_JOBS[$I_PARTICLE]}`; do
		    FILE_ID="Job"$JOB_ID"_"`date +%Y%m%d`"_"`date +%H%M%S`
		    
		    logfile=$OUTDIR/"output_"$FILE_ID".log"
		    
		    new_fclfile=$OUTDIR"/ProtonAbsorberRuns_"${PARTICLES[$I_PARTICLE]}"_"$FILE_ID".fcl"
		    cat ~/Mu2e/Offline/Projects/edmonds/IPASims/fcl/ProtonAbsorberRuns_${PARTICLES[$I_PARTICLE]}.fcl > $new_fclfile
		    seed=$RANDOM
		    echo -e "\nservices.user.SeedService.baseSeed         :  "$seed >> $new_fclfile
		    
		    new_rootfile=$OUTDIR"/PAbsRun_"${PARTICLES[$I_PARTICLE]}"_"$FILE_ID".root"
		    echo -e "\nservices.TFileService.fileName  : \""$new_rootfile"\"" >> $new_fclfile

		    let TOTAL_JOBS=$TOTAL_JOBS+1;
		    echo $new_rootfile

		    nohup mu2e --config $new_fclfile --nevts ${N_EVTS[$I_PARTICLE]} > $logfile 2>&1 &
		done # jobs per particle
	    done # particles

	    while true ; do
		TEST=`ps -u aedmonds | grep mu2e | wc -l`
		if [ $TEST -eq 0 ] ; then
		    break;
		else
#		    echo "Waiting..."
		    sleep 10s
		fi
	    done
	done
    done
done

echo "Done!"