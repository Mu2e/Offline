i=0
layerOffset=42
moduleGap=5
#  for moduleGap in {1..5}
#  do

    dz=`echo "804.3+$moduleGap" | bc -l`

#    for layerOffset in {0..62..2}
#    do
#      for photonYield in {3189,3781,4359}  # 24,28,32 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
      for photonYield in {2604,2897,3189,3481,3774,4066,4359,4651,4943,5236,5528,5820,6113,6405,6697,6990,7282,7574,7867}  # 20,22,...,56 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
      do

        ((i++));

        overlayFactor="0"

        genconfigfile=CRVResponse/efficiencyCheck/submit/genconfig_5cm_6600'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/genconfig_5cm_6600.txt\"" >| $genconfigfile
        echo "double cosmicDYB.dz = $dz;" >> $genconfigfile

        geomfile=CRVResponse/efficiencyCheck/submit/geom_5cm_6600'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/geom_5cm_6600.txt\"" >| $geomfile
        echo "double crs.gapBetweenModules = $moduleGap;" >> $geomfile
        echo "double crs.layerOffset       = $layerOffset;" >> $geomfile

#### for adjacentPulseTimeDifference = 5ns

#        name=CRV_Efficiency_check_5cm_6600
#        fclfile=CRVResponse/efficiencyCheck/submit/$name'_'$i.fcl
#        fcllist=CRVResponse/efficiencyCheck/submit/$name'_'$i.txt
#        echo "#include \"CRVResponse/efficiencyCheck/CRV_Efficiency_check_5cm_6600.fcl\"" >| $fclfile
#        echo "services.GeometryService.inputFile                      : \"$geomfile\"" >> $fclfile
#        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
#        echo "physics.producers.CrvPhotonArrivals.scintillationYield  : $photonYield" >> $fclfile
#        echo "physics.producers.backgroundOverlay.overlayFactor       : $overlayFactor" >> $fclfile

#        generate_fcl --description=$name --dsconf=$i --run=1 --events=20000 --njobs=50 $fclfile
#        ls $PWD/000/cnf.$USER.$name.$i.*.fcl > $fcllist
#        clustername=$name'_'gap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield
#        mu2eprodsys --setup=./setup.sh --fcllist=$fcllist --clustername=$clustername --dsconf=$i --wfproject=$name

#### for adjacentPulseTimeDifference = 10ns

        name=CRV_Efficiency_check_5cm10_6600
        fclfile=CRVResponse/efficiencyCheck/submit/$name'_'$i.fcl
        fcllist=CRVResponse/efficiencyCheck/submit/$name'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/CRV_Efficiency_check_5cm10_6600.fcl\"" >| $fclfile
        echo "services.GeometryService.inputFile                      : \"$geomfile\"" >> $fclfile
        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
        echo "physics.producers.CrvPhotonArrivals.scintillationYield  : $photonYield" >> $fclfile
        echo "physics.producers.backgroundOverlay.overlayFactor       : $overlayFactor" >> $fclfile

        generate_fcl --description=$name --dsconf=$i --run=1 --events=20000 --njobs=50 $fclfile
        ls $PWD/000/cnf.$USER.$name.$i.*.fcl > $fcllist
        clustername=$name'_'gap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield
        mu2eprodsys --setup=./setup.sh --fcllist=$fcllist --clustername=$clustername --dsconf=$i --wfproject=$name

      done
#    done
#  done
