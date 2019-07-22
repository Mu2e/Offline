#To submit the jobs:
#setup gridexport
#Delete content of CRVResponse/efficiency/submit
#Delete all lines of CRVResponse/efficiencyCheck/jobs.sh except first line
#Run this script from workspace/Offline
#Run gridexport from workspace/Offline and copy the location of the zip file to the joblist script to $CODELOCATION
#Run the joblist script from workspace/Offline

sed '378s/.*/CrvCoincidenceX  : @local::CrvCoincidenceT/' CRVResponse/efficiencyCheck/CRV_Efficiency_check.fcl >| CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_T.fcl
sed '378s/.*/CrvCoincidenceX  : @local::CrvCoincidenceTS/' CRVResponse/efficiencyCheck/CRV_Efficiency_check.fcl >| CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_TS.fcl
sed '378s/.*/CrvCoincidenceX  : @local::CrvCoincidenceE/' CRVResponse/efficiencyCheck/CRV_Efficiency_check.fcl >| CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_E.fcl
sed '378s/.*/CrvCoincidenceX  : @local::CrvCoincidenceR/' CRVResponse/efficiencyCheck/CRV_Efficiency_check.fcl >| CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_R.fcl
sed '378s/.*/CrvCoincidenceX  : @local::CrvCoincidenceL/' CRVResponse/efficiencyCheck/CRV_Efficiency_check.fcl >| CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_L.fcl
sed '378s/.*/CrvCoincidenceX  : @local::CrvCoincidenceU/' CRVResponse/efficiencyCheck/CRV_Efficiency_check.fcl >| CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_U.fcl
sed '378s/.*/CrvCoincidenceX  : @local::CrvCoincidenceD/' CRVResponse/efficiencyCheck/CRV_Efficiency_check.fcl >| CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_D.fcl

sed -i '440s/.*/         pars: [@local::outsideCrvCutT, @local::belowCrvCutT]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_T.fcl
sed -i '440s/.*/         pars: [@local::outsideCrvCutT, @local::belowCrvCutT]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_TS.fcl
sed -i '440s/.*/         pars: [@local::outsideCrvCutE, @local::belowCrvCutE]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_E.fcl
sed -i '440s/.*/         pars: [@local::outsideCrvCutR, @local::belowCrvCutR]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_R.fcl
sed -i '440s/.*/         pars: [@local::outsideCrvCutL, @local::belowCrvCutL]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_L.fcl
sed -i '440s/.*/         pars: [@local::outsideCrvCutU, @local::belowCrvCutU]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_U.fcl
sed -i '440s/.*/         pars: [@local::outsideCrvCutD, @local::belowCrvCutD]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_D.fcl

sed -i '690s/.*/physics.producers.CrvPhotons.lookupTableFileNames : ["CRVConditions\/v6_0\/LookupTable_6000_0"]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_T.fcl
sed -i '690s/.*/physics.producers.CrvPhotons.lookupTableFileNames : ["CRVConditions\/v6_0\/LookupTable_6000_1"]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_TS.fcl
sed -i '690s/.*/physics.producers.CrvPhotons.lookupTableFileNames : ["CRVConditions\/v6_0\/LookupTable_5000_1"]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_E.fcl
sed -i '690s/.*/physics.producers.CrvPhotons.lookupTableFileNames : ["CRVConditions\/v6_0\/LookupTable_4550_0"]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_R.fcl
sed -i '690s/.*/physics.producers.CrvPhotons.lookupTableFileNames : ["CRVConditions\/v6_0\/LookupTable_4550_0"]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_L.fcl
sed -i '690s/.*/physics.producers.CrvPhotons.lookupTableFileNames : ["CRVConditions\/v6_0\/LookupTable_6900_1"]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_U.fcl
sed -i '690s/.*/physics.producers.CrvPhotons.lookupTableFileNames : ["CRVConditions\/v6_0\/LookupTable_5700_0"]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_D.fcl

sed -i '691s/.*/physics.producers.CrvPhotons.reflectors           : [ 0 ]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_T.fcl
sed -i '691s/.*/physics.producers.CrvPhotons.reflectors           : [ 1 ]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_TS.fcl
sed -i '691s/.*/physics.producers.CrvPhotons.reflectors           : [-1 ]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_E.fcl
sed -i '691s/.*/physics.producers.CrvPhotons.reflectors           : [ 0 ]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_R.fcl
sed -i '691s/.*/physics.producers.CrvPhotons.reflectors           : [ 0 ]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_L.fcl
sed -i '691s/.*/physics.producers.CrvPhotons.reflectors           : [ 1 ]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_U.fcl
sed -i '691s/.*/physics.producers.CrvPhotons.reflectors           : [ 0 ]/' CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_D.fcl

i=0
name=CRV_Efficiency_check_L
layerOffset=42
#layerOffset=0
moduleGap="3.0"
gapSmall="0.0"
gapLarge="0.46"
#PEYield=44
joblist=CRVResponse/efficiencyCheck/jobs.sh

#for layerOffset in {0..62..1}
#do
    for PEYield in {10..68..2}  #PE per SiPM @ 1m away from SiPM of 3m long counter
    do

        ((i++));

        photonYield=`echo "39400.0*$PEYield/68.0" | bc -l`
        dx=`echo "821.44+$moduleGap+7*$gapLarge+8*$gapSmall" | bc -l`
        dy=`echo "821.44+$moduleGap+7*$gapLarge+8*$gapSmall" | bc -l`
        dz=`echo "821.44+$moduleGap+7*$gapLarge+8*$gapSmall" | bc -l`

        genconfigfile=CRVResponse/efficiencyCheck/submit/genconfig'_'$name'_'$i.txt
        #echo "#include \"CRVResponse/efficiencyCheck/genconfig_T.txt\"" >| $genconfigfile
        #echo "#include \"CRVResponse/efficiencyCheck/genconfig_T4.txt\"" >| $genconfigfile
        #echo "#include \"CRVResponse/efficiencyCheck/genconfig_TS.txt\"" >| $genconfigfile
        #echo "#include \"CRVResponse/efficiencyCheck/genconfig_E.txt\"" >| $genconfigfile
        #echo "#include \"CRVResponse/efficiencyCheck/genconfig_R.txt\"" >| $genconfigfile
        echo "#include \"CRVResponse/efficiencyCheck/genconfig_L.txt\"" >| $genconfigfile
        #echo "#include \"CRVResponse/efficiencyCheck/genconfig_U.txt\"" >| $genconfigfile
        #echo "#include \"CRVResponse/efficiencyCheck/genconfig_D.txt\"" >| $genconfigfile

        #echo "double cosmicFromTH2.dx = $dx;" >> $genconfigfile  #for E
        #echo "double cosmicFromTH2.dy = $dy;" >> $genconfigfile  #for U,D
        echo "double cosmicFromTH2.dz = $dz;" >> $genconfigfile  #for T,T4,TS,R,L

        geomfile=CRVResponse/efficiencyCheck/submit/geom'_'$name'_'$i.txt
        #echo "#include \"CRVResponse/efficiencyCheck/geom_T.txt\"" >| $geomfile
        #echo "#include \"CRVResponse/efficiencyCheck/geom_TS.txt\"" >| $geomfile
        #echo "#include \"CRVResponse/efficiencyCheck/geom_E.txt\"" >| $geomfile
        #echo "#include \"CRVResponse/efficiencyCheck/geom_R.txt\"" >| $geomfile
        echo "#include \"CRVResponse/efficiencyCheck/geom_L.txt\"" >| $geomfile
        #echo "#include \"CRVResponse/efficiencyCheck/geom_U.txt\"" >| $geomfile
        #echo "#include \"CRVResponse/efficiencyCheck/geom_D.txt\"" >| $geomfile

        echo "double crs.gapBetweenModules = $moduleGap;" >> $geomfile
        echo "double crs.gapSmall          = $gapSmall;" >> $geomfile
        echo "double crs.gapLarge          = $gapLarge;" >> $geomfile
        echo "double crs.layerOffset       = $layerOffset;" >> $geomfile

        fclfile=CRVResponse/efficiencyCheck/submit/$name'_'$i.fcl
        #echo "#include \"CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_T.fcl\"" >| $fclfile
        #echo "#include \"CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_TS.fcl\"" >| $fclfile
        #echo "#include \"CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_E.fcl\"" >| $fclfile
        #echo "#include \"CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_R.fcl\"" >| $fclfile
        echo "#include \"CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_L.fcl\"" >| $fclfile
        #echo "#include \"CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_U.fcl\"" >| $fclfile
        #echo "#include \"CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_D.fcl\"" >| $fclfile

        echo "services.GeometryService.inputFile                      : \"$geomfile\"" >> $fclfile
        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
        echo "physics.producers.CrvPhotons.scintillationYield         : $photonYield" >> $fclfile
        echo "physics.producers.CrvSiPMCharges.deadSiPMProbability    : 0.02" >> $fclfile

        generate_fcl --description=$name --dsconf=$i --run=1 --events=100000 --njobs=20 --include $fclfile
        tar -zcvf fcls.$name.$i.tar.gz 000/cnf.$USER.$name.$i.*.fcl
        rm -rf 000
#        mv fcls.$name.$i.tar.gz /pnfs/mu2e/scratch/outstage/ehrlich/fcls/.
        clustername=$name'_'moduleGap$moduleGap'_'gapLarge$gapLarge'_'gapSmall$gapSmall'_'layerOffset$layerOffset'_'PEYieldYield$PEYield
        echo "mu2eprodsys --expected-lifetime=8h --code=\$CODELOCATION --fcllist=/pnfs/mu2e/scratch/outstage/ehrlich/fcls/fcls.$name.$i.tar.gz --clustername=$clustername --dsconf=$i --wfproject=$name" >> $joblist

    done
#done
