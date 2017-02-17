layerOffset=42
timeWindow=10
sides=2
PEthreshold=12
moduleGap=3

    for photonYield in {3156,3738,4320}  # 24,28,32 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
    do
      for layerOffset in {0..62..2}
      do

              directory=/pnfs/mu2e/scratch/users/ehrlich/workflow
              files=`ls $directory/CRV_Efficiency_check_5cm/outstage/*.CRV_efficiency5cm_top_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield/*/*/log.*.log`
#              files=`ls $directory/CRV_Efficiency_check_5cm10/outstage/*.CRV_efficiency5cm10_top_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield/*/*/log.*.log`
#              files=`ls $directory/CRV_Efficiency_check_5cm_6600/outstage/*.CRV_efficiency5cm_top6600_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield/*/*/log.*.log`
#              files=`ls $directory/CRV_Efficiency_check_5cm10_6600/outstage/*.CRV_efficiency5cm10_top6600_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield/*/*/log.*.log`
#              files=`ls $directory/CRV_Efficiency_check_5cm_verticalPlanes/outstage/*.CRV_efficiency5cm_side_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield/*/*/log.*.log`
#              files=`ls $directory/CRV_Efficiency_check_5cm_downstreamPlanes/outstage/*.CRV_efficiency5cm_downstream_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield/*/*/log.*.log`

              events=0
              eventsCoincidence=0

              for file in $files;
              do

                searchString="SUMMARY CrvCoincidencePE$PEthreshold""T$timeWindow"
                searchCommand="grep '$searchString' $file"
#                searchCommand="tail -n 70 $file | grep '$searchString'"
                searchResult=`eval $searchCommand`
                if [ $? -ne 0 ]; then
                  continue
                fi

                s=($searchResult)
                eventsTmp=${s[4]}
                eventsCoincidenceTmp=${s[2]}

                events=$((events+eventsTmp))
                eventsCoincidence=$((eventsCoincidence+eventsCoincidenceTmp))

              done

              efficiency=`echo "$eventsCoincidence/$events" | bc -l`
              echo modulegap:$moduleGap layerOffset:$layerOffset photonYield:$photonYield PEThreshold:$PEthreshold timeWindow:$timeWindow sidesChecked:$sides eventsCoincidence:$eventsCoincidence eventsTotal:$events efficiency:$efficiency

      done
    done
