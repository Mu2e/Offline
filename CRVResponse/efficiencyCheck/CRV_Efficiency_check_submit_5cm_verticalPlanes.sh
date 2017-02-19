i=0
layerOffset=42
moduleGap=3
#  for moduleGap in {1..5}
#  do

    dz=$((811+$moduleGap))

#    for layerOffset in {0..62..2}
#    do
#      for photonYield in {3156,3738,4320}  # 24,28,32 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
      for photonYield in {2575,2865,3156,3447,3738,4029,4320,4611,4901,5192,5483,5774,6065,6356,6646,6937,7228,7519,7810}  # 20,22,...,56 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
      do

        ((i++));

        overlayFactor="0"

        genconfigfile=CRVResponse/efficiencyCheck/submit/genconfig_5cm_verticalPlanes'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/genconfig_5cm_verticalPlanes.txt\"" >| $genconfigfile
        echo "double cosmicDYB.dz = $dz;" >> $genconfigfile

        geomfile=CRVResponse/efficiencyCheck/submit/geom_5cm_verticalPlanes'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/geom_5cm_verticalPlanes.txt\"" >| $geomfile
        echo "double crs.gapBetweenModules = $moduleGap;" >> $geomfile
        echo "double crs.layerOffset       = $layerOffset;" >> $geomfile

        name=CRV_Efficiency_check_5cm_verticalPlanes
        fclfile=CRVResponse/efficiencyCheck/submit/$name'_'$i.fcl
        fcllist=CRVResponse/efficiencyCheck/submit/$name'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/CRV_Efficiency_check_5cm_verticalPlanes.fcl\"" >| $fclfile
        echo "services.user.GeometryService.inputFile                 : \"$geomfile\"" >> $fclfile
        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
        echo "physics.producers.CrvPhotonArrivals.scintillationYield  : $photonYield" >> $fclfile
        echo "physics.producers.backgroundOverlay.overlayFactor       : $overlayFactor" >> $fclfile

        generate_fcl --description=$name --dsconf=$i --run=1 --events=20000 --njobs=50 $fclfile
        ls /mu2e/app/users/ehrlich/work_08302015/Offline/000/cnf.ehrlich.$name.$i.*.fcl > $fcllist
        clustername=CRV_efficiency5cm_side_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield
        mu2eprodsys --setup=./setup.sh --fcllist=$fcllist --clustername=$clustername --dsconf=$i --wfproject=$name

      done
#    done
#  done
