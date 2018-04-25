name=CRV_Efficiency_check_5cm_6000_e
layerOffset=42
moduleGap="3.0"
gapSmall="0.1"
gapLarge="2.0"

timeWindow=10
sides=2
#PEthreshold=12

for PEthreshold in {10..40..2}
do
    for PEYield in {32..68..2}  #PE per SiPM @ 1m away from SiPM of 3m long counter
    do

              directory=/pnfs/mu2e/scratch/users/ehrlich/workflow

              files=$directory/$name/outstage/'*.'$name'_'moduleGap$moduleGap'_'gapLarge$gapLarge'_'gapSmall$gapSmall'_'layerOffset$layerOffset'_'PEYieldYield$PEYield/*/*/log.*.log

              events=0
              eventsCoincidence=0

              for file in $files;
              do

                searchString="SUMMARY CrvCoincidencePE$PEthreshold""T$timeWindow"
#                searchCommand="grep '$searchString' $file"
                searchCommand="timeout 10 tail -n 200 $file | grep '$searchString'"
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
              echo modulegap:$moduleGap layerOffset:$layerOffset PEYield:$PEYield PEThreshold:$PEthreshold timeWindow:$timeWindow sidesChecked:$sides eventsCoincidence:$eventsCoincidence eventsTotal:$events efficiency:$efficiency

    done
done
