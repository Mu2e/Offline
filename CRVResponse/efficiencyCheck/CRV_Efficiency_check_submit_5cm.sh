i=0
layerOffset=42
moduleGap=3
#  for moduleGap in {1..5}
#  do

    dz=$((811+$moduleGap))

#    for layerOffset in {0..62..2}
#    do
#      for photonYield in {3156,3738,4320}  # 24,28,32 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
      for photonYield in {2575,2865,3156,3447,3738, 4029,4320,4611,4901,5192, 5483,5774,6065,6356,6646, 6937,7228,7519,7810}  # 20,22,...,56 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
      do

        ((i++));

        overlayFactor="0"

        genconfigfile=CRVResponse/efficiencyCheck/submit/genconfig_5cm'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/genconfig_5cm.txt\"" >| $genconfigfile
        echo "double cosmicDYB.dz = $dz;" >> $genconfigfile

        geomfile=CRVResponse/efficiencyCheck/submit/geom_5cm'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/geom_5cm.txt\"" >| $geomfile
        echo "double crs.gapBetweenModules = $moduleGap;" >> $geomfile
        echo "double crs.layerOffset       = $layerOffset;" >> $geomfile

#### for adjacentPulseTimeDifference = 5ns

        fclfile=CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_5cm'_'$i.fcl
        echo "#include \"CRVResponse/efficiencyCheck/CRV_Efficiency_check_5cm.fcl\"" >| $fclfile
        echo "services.user.GeometryService.inputFile                 : \"$geomfile\"" >> $fclfile
        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
        echo "physics.producers.CrvPhotonArrivals.scintillationYield  : $photonYield" >> $fclfile
        echo "physics.producers.backgroundOverlay.overlayFactor       : $overlayFactor" >> $fclfile

        mu2eart --setup=./setup.sh --fcl=$fclfile --njobs=100 --events-per-job=10000 --jobname=CRV_efficiency5cm_top_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield --outstage=/pnfs/mu2e/scratch/outstage

#### for adjacentPulseTimeDifference = 10ns

        fclfile=CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_5cm10'_'$i.fcl
        echo "#include \"CRVResponse/efficiencyCheck/CRV_Efficiency_check_5cm10.fcl\"" >| $fclfile
        echo "services.user.GeometryService.inputFile                 : \"$geomfile\"" >> $fclfile
        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
        echo "physics.producers.CrvPhotonArrivals.scintillationYield  : $photonYield" >> $fclfile
        echo "physics.producers.backgroundOverlay.overlayFactor       : $overlayFactor" >> $fclfile

#        mu2eart --setup=./setup.sh --fcl=$fclfile --njobs=50 --events-per-job=20000 --jobname=CRV_efficiency5cm10_top_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield --outstage=/pnfs/mu2e/scratch/outstage

      done
#    done
#  done
