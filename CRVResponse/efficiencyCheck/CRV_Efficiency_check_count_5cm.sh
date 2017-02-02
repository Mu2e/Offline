layerOffset=42
timeWindow=10
sides=2
PEthreshold=10
moduleGap=3

    for photonYield in {3047,3555,4063,4570}  # 24,28,32,36 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
    do
      for layerOffset in {0..62..2}
      do

              files=`ls /pnfs/mu2e/scratch/outstage/ehrlich/CRV_efficiency5cm_top_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield.*/*/mu2e.log`
#              files=`ls /pnfs/mu2e/scratch/outstage/ehrlich/CRV_efficiency5cm10_top_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield.*/*/mu2e.log`
#              files=`ls /pnfs/mu2e/scratch/outstage/ehrlich/CRV_efficiency5cm_top6600_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield.*/*/mu2e.log`
#              files=`ls /pnfs/mu2e/scratch/outstage/ehrlich/CRV_efficiency5cm10_top6600_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield.*/*/mu2e.log`
#              files=`ls /pnfs/mu2e/scratch/outstage/ehrlich/CRV_efficiency5cm_side_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield.*/*/mu2e.log`
#              files=`ls /pnfs/mu2e/scratch/outstage/ehrlich/CRV_efficiency5cm_downstream_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield.*/*/mu2e.log`

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
