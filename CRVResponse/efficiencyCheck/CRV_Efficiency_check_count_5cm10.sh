layerOffset=42
moduleGap=5
timeWindow=10
sides=2
#  for moduleGap in {2..5}
#  do
#    for layerOffset in {0..62..2}
#    do

       for PEthreshold in {6..16..2}
       do

         for photonYield in {3000,3500,4000,4500,5000,5500,6000,6500,7000}  # 20,??,26,??,33,??,40,??,?? PE/SiPM for 5cm wide / 5.6m long counter
         do

              files=`ls /pnfs/mu2e/scratch/outstage/ehrlich/CRV_efficiency5cm10_top_moduleGap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield.*/*/mu2e.log`

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

              inefficiency=`echo "$eventsCoincidence/$events" | bc -l`
              echo modulegap:$moduleGap layerOffset:$layerOffset photonYield:$photonYield PEThreshold:$PEthreshold timeWindow:$timeWindow sidesChecked:$sides eventsCoincidence:$eventsCoincidence eventsTotal:$events inefficiency:$inefficiency

         done
       done
#    done
#  done
