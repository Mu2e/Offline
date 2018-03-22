i=0
layerOffset=0
moduleGap="4.5"
joblist=CRVResponse/efficiencyCheck/jobs2.sh
#  for moduleGap in {"4.5","6.5"}
#  do

    dz=`echo "804.3+$moduleGap" | bc -l`
echo "*****************************"
echo "WIDER INNER AND OUTER GAP!!!!"
echo "*****************************"
    dz=`echo "818.5+$moduleGap" | bc -l`

#    for layerOffset in {0..62..2}
#    do
#      for photonYield in {3189,3781,4359}  # 24,28,32 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
      for photonYield in {2604,2897,3189,3481,3774,4066,4359,4651,4943,5236,5528,5820,6113,6405,6697,6990,7282,7574,7867}  # 20,22,...,56 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
      do

        ((i++));

        overlayFactor="0"

        genconfigfile=CRVResponse/efficiencyCheck/submit/genconfig_5cm_5000r'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/genconfig_5cm_5000.txt\"" >| $genconfigfile
        echo "double cosmicDYB.dz = $dz;" >> $genconfigfile

        geomfile=CRVResponse/efficiencyCheck/submit/geom_5cm_5000r'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/geom_5cm_5000.txt\"" >| $geomfile
        echo "double crs.gapBetweenModules = $moduleGap;" >> $geomfile
        echo "double crs.gapSmall          = 0.5;" >> $geomfile
        echo "double crs.gapLarge          = 1.0;" >> $geomfile
        echo "double crs.layerOffset       = $layerOffset;" >> $geomfile

#### for adjacentPulseTimeDifference = 0ns

#        name=CRV_Efficiency_check_5cm0_5000r
#        fclfile=CRVResponse/efficiencyCheck/submit/$name'_'$i.fcl
#        fcllist=CRVResponse/efficiencyCheck/submit/$name'_'$i.txt
#        echo "#include \"CRVResponse/efficiencyCheck/CRV_Efficiency_check_5cm0_5000r.fcl\"" >| $fclfile
#        echo "services.GeometryService.inputFile                      : \"$geomfile\"" >> $fclfile
#        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
#        echo "physics.producers.CrvPhotons.scintillationYield         : $photonYield" >> $fclfile
#        echo "physics.producers.backgroundOverlay.overlayFactor       : $overlayFactor" >> $fclfile
#
#        generate_fcl --description=$name --dsconf=$i --run=1 --events=20000 --njobs=50 $fclfile
#        ls $PWD/000/cnf.$USER.$name.$i.*.fcl > $fcllist
#        clustername=$name'_'gap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield
#        mu2eprodsys --setup=./setup.sh --fcllist=$fcllist --clustername=$clustername --dsconf=$i --wfproject=$name

#### for adjacentPulseTimeDifference = 5ns

        name=CRV_Efficiency_check_5cm_5000r
        fclfile=CRVResponse/efficiencyCheck/submit/$name'_'$i.fcl
        echo "#include \"CRVResponse/efficiencyCheck/CRV_Efficiency_check_5cm_5000r.fcl\"" >| $fclfile
        echo "services.GeometryService.inputFile                      : \"$geomfile\"" >> $fclfile
        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
        echo "physics.producers.CrvPhotons.scintillationYield         : $photonYield" >> $fclfile
        echo "physics.producers.backgroundOverlay.overlayFactor       : $overlayFactor" >> $fclfile

        generate_fcl --description=$name --dsconf=$i --run=1 --events=50000 --njobs=20 $fclfile
        tar -zcvf fcls.$name.$i.tar.gz 000/cnf.$USER.$name.$i.*.fcl
        rm -rf 000
        mv fcls.$name.$i.tar.gz /pnfs/mu2e/scratch/outstage/ehrlich/fcls/.
        clustername=$name'_'gap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield
        echo "mu2eprodsys --code=\$CODELOCATION --fcllist=/pnfs/mu2e/scratch/outstage/ehrlich/fcls/fcls.$name.$i.tar.gz --clustername=$clustername --dsconf=$i --wfproject=$name" >> $joblist

#### for adjacentPulseTimeDifference = 10ns

#        name=CRV_Efficiency_check_5cm10_5000r
#        fclfile=CRVResponse/efficiencyCheck/submit/$name'_'$i.fcl
#        fcllist=CRVResponse/efficiencyCheck/submit/$name'_'$i.txt
#        echo "#include \"CRVResponse/efficiencyCheck/CRV_Efficiency_check_5cm10_5000r.fcl\"" >| $fclfile
#        echo "services.GeometryService.inputFile                      : \"$geomfile\"" >> $fclfile
#        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
#        echo "physics.producers.CrvPhotons.scintillationYield         : $photonYield" >> $fclfile
#        echo "physics.producers.backgroundOverlay.overlayFactor       : $overlayFactor" >> $fclfile
#
#        generate_fcl --description=$name --dsconf=$i --run=1 --events=20000 --njobs=50 $fclfile
#        ls $PWD/000/cnf.$USER.$name.$i.*.fcl > $fcllist
#        clustername=$name'_'gap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield
#        mu2eprodsys --setup=./setup.sh --fcllist=$fcllist --clustername=$clustername --dsconf=$i --wfproject=$name

      done
#    done
#  done
