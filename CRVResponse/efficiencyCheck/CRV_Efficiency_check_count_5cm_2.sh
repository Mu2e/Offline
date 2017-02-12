layerOffset=42
timeWindow=10
sides=2
PEthreshold=12
moduleGap=3

for PEthreshold in {10..20..2}
do
    for photonYield in {3047,3301,3555,3809,4063,4316,4570,4824,5078,5332,5586,5840,6094,6348}  # 24,26,28,30,32,34,36,38,40,42,44,46,48,50 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
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
