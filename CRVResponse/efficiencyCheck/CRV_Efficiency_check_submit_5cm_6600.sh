i=0
layerOffset=42
moduleGap=3
#  for moduleGap in {1..5}
#  do

    dz=$((811+$moduleGap))

    for layerOffset in {0..62..2}
    do
      for photonYield in {3047,3555,4063,4570}  # 24,28,32,36 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
      do

        ((i++));

        scintTolerance="0.2"

        genconfigfile=CRVResponse/efficiencyCheck/submit/genconfig_5cm_6600'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/genconfig_5cm_6600.txt\"" >| $genconfigfile
        echo "double cosmicDYB.dz = $dz;" >> $genconfigfile

        geomfile=CRVResponse/efficiencyCheck/submit/geom_5cm_6600'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/geom_5cm_6600.txt\"" >| $geomfile
        echo "double crs.gapBetweenModules = $moduleGap;" >> $geomfile
        echo "double crs.layerOffset       = $layerOffset;" >> $geomfile

#### for adjacentPulseTimeDifference = 5ns

        fclfile=CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_5cm_6600'_'$i.fcl
        echo "#include \"CRVResponse/efficiencyCheck/CRV_Efficiency_check_5cm_6600.fcl\"" >| $fclfile
        echo "services.user.GeometryService.inputFile                 : \"$geomfile\"" >> $fclfile
        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
        echo "physics.producers.CrvPhotonArrivals.scintillationYield  : $photonYield" >> $fclfile
        echo "physics.producers.CrvPhotonArrivals.scintillationYieldTolerance  : $scintTolerance" >> $fclfile
#        echo "physics.producers.CrvSiPMResponses.ThermalProb          : 0" >> $fclfile

        mu2eart --setup=./setup.sh --fcl=$fclfile --njobs=50 --events-per-job=20000 --jobname=CRV_efficiency5cm_top6600_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield --outstage=/pnfs/mu2e/scratch/outstage

#### for adjacentPulseTimeDifference = 10ns

        fclfile=CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_5cm10_6600'_'$i.fcl
        echo "#include \"CRVResponse/efficiencyCheck/CRV_Efficiency_check_5cm10_6600.fcl\"" >| $fclfile
        echo "services.user.GeometryService.inputFile                 : \"$geomfile\"" >> $fclfile
        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
        echo "physics.producers.CrvPhotonArrivals.scintillationYield  : $photonYield" >> $fclfile
        echo "physics.producers.CrvPhotonArrivals.scintillationYieldTolerance  : $scintTolerance" >> $fclfile
#        echo "physics.producers.CrvSiPMResponses.ThermalProb          : 0" >> $fclfile

        mu2eart --setup=./setup.sh --fcl=$fclfile --njobs=50 --events-per-job=20000 --jobname=CRV_efficiency5cm10_top6600_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield --outstage=/pnfs/mu2e/scratch/outstage

      done
    done
#  done
