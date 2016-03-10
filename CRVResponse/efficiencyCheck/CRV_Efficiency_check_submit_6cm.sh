i=0
layerOffset=42
moduleGap=5
#  for moduleGap in {2..5}
#  do

     dz=$(echo 727.72 + $moduleGap | bc)

#    for layerOffset in {0..62..2}
#    do
      for photonYield in {3000,3500,4000,4500,5000,5500,6000,6500,7000}  # 17,??,22,??,28,??,33,??,?? PE/SiPM for 6cm wide / 5.6m long counter
      do

        ((i++));

        scintTolerance="0"
#        backgroundParam1="2.9e-4"
        backgroundParam1="0"
        backgroundParam3=$(echo "scale=3; 5000/$photonYield" | bc)

        genconfigfile=CRVResponse/efficiencyCheck/submit/genconfig_6cm'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/genconfig_6cm.txt\"" >| $genconfigfile
        echo "double cosmicDYB.dz = $dz;" >> $genconfigfile

        geomfile=CRVResponse/efficiencyCheck/submit/geom_6cm'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/geom_6cm.txt\"" >| $geomfile
        echo "double crs.gapBetweenModules = $moduleGap;" >> $geomfile
        echo "double crs.layerOffset       = $layerOffset;" >> $geomfile

        fclfile=CRVResponse/efficiencyCheck/submit/CRV_Efficiency_check_6cm'_'$i.fcl
        echo "#include \"CRVResponse/efficiencyCheck/CRV_Efficiency_check_6cm.fcl\"" >| $fclfile
        echo "services.user.GeometryService.inputFile                 : \"$geomfile\"" >> $fclfile
        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
        echo "physics.producers.CrvPhotonArrivals.scintillationYield  : $photonYield" >> $fclfile
        echo "physics.producers.CrvPhotonArrivals.scintillationYieldTolerance  : $scintTolerance" >> $fclfile
        echo "physics.producers.CrvSiPMResponses.ThermalProb          : 0" >> $fclfile
        echo "physics.producers.CrvSiPMResponses.backgroundParam1     : $backgroundParam1" >> $fclfile
        echo "physics.producers.CrvSiPMResponses.backgroundParam3     : $backgroundParam3" >> $fclfile

        mu2eart --setup=./setup.sh --fcl=$fclfile --njobs=50 --events-per-job=10000 --jobname=CRV_efficiency6cm_top_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield --outstage=/pnfs/mu2e/scratch/outstage

      done
#    done
#  done
