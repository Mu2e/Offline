#To submit the jobs:
setup gridexport
setup mu2etools
#Delete content of CRVResponse/efficiency/submit
#Delete all lines of CRVResponse/efficiencyCheck/jobs.sh except first line
#Run this script from workspace/Offline
#Run gridexport from workspace/Offline and copy the location of the zip file to the joblist script to $CODELOCATION
#Run the joblist script from workspace/Offline

sed '378s/.*/CrvCoincidenceX  : @local::CrvCoincidenceL/' Mu2eANL/Analyses/fcl/CRV_Efficiency_check.fcl >| CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_L_teststand.fcl
#sed -i '441s/.*/         pars: [@local::outsideCrvCutT, @local::belowCrvCutT]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_L_teststand.fcl
#sed -i '691s/.*/physics.producers.CrvPhotons.lookupTableFileNames : ["CRVConditions\/v6_0\/LookupTable_4550_0"]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_L_teststand.fcl
#sed -i '692s/.*/physics.producers.CrvPhotons.reflectors           : [ 0 ]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_L_teststand.fcl

i=0
name=CRV_Efficiency_check_L_teststand0
layerOffset=42
#layerOffset=0
moduleGap="3.0"
gapSmall="0.0"
gapLarge="0.46"
#PEYield=44
joblist=CRVResponse/efficiencyCheck/jobs.sh

#for layerOffset in {0..62..1}
#do
    for PEYield in {40..42..2}  #PE per SiPM @ 1m away from SiPM of 3m long counter - Originally {10..68..2}
    do
        ((i++));
        photonYield=`echo "39400.0*$PEYield/68.0" | bc -l`
        dx=`echo "821.44+$moduleGap+7*$gapLarge+8*$gapSmall" | bc -l`
        dy=`echo "821.44+$moduleGap+7*$gapLarge+8*$gapSmall" | bc -l`
        dz=`echo "821.44+$moduleGap+7*$gapLarge+8*$gapSmall" | bc -l`

        genconfigfile=CRVResponse/efficiencyCheck/submit/genconfig'_'$name'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/genconfig_L_teststand.txt\"" >| $genconfigfile
        echo "double cosmicDYB.dz = $dz;" >> $genconfigfile  #for teststand for T,T4,TS,R,L
        geomfile=CRVResponse/efficiencyCheck/submit/geom'_'$name'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/geom_L_teststand0.txt\"" >| $geomfile
        echo "double crs.gapBetweenModules = $moduleGap;" >> $geomfile
        echo "double crs.gapSmall          = $gapSmall;" >> $geomfile
        echo "double crs.gapLarge          = $gapLarge;" >> $geomfile
        echo "double crs.layerOffset       = $layerOffset;" >> $geomfile

        fclfile=CRVResponse/efficiencyCheck/submit/$name'_'$i.fcl
        echo "#include \"CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_L_teststand.fcl\"" >| $fclfile
        echo "services.GeometryService.inputFile                      : \"$geomfile\"" >> $fclfile
        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
        echo "physics.producers.CrvPhotons.scintillationYield         : $photonYield" >> $fclfile
        echo "physics.producers.CrvSiPMCharges.deadSiPMProbability    : 0.00" >> $fclfile

	#Originally 100000 events
        generate_fcl --description=$name --dsconf=$i --run=1 --events=100000 --njobs=100 --include $fclfile
        tar -zcvf fcls.$name.$i.tar.gz 000/cnf.$USER.$name.$i.*.fcl
        rm -rf 000
        DATE=$(date +"%Y%m%d%H%M%S")
        OUTDIR=/pnfs/mu2e/scratch/users/${USER}/fcl/${DATE}/
        mkdir -p ${OUTDIR}
        mv fcls.$name.$i.tar.gz ${OUTDIR}
        clustername=$name'_'moduleGap$moduleGap'_'gapLarge$gapLarge'_'gapSmall$gapSmall'_'layerOffset$layerOffset'_'PEYieldYield$PEYield
        echo "mu2eprodsys --expected-lifetime=8h --memory=4GB --code=\$CODELOCATION --fcllist=${OUTDIR}/fcls.$name.$i.tar.gz --clustername=$clustername --dsconf=$i --wfproject=$name" >> $joblist

    done
#done
