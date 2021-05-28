name=CRV_Efficiency_check_T
layerOffset=42
#layerOffset=0
moduleGap="2.0"
gapSmall="0.0"
gapLarge="0.46"

timeWindow=10
#timeWindow=20
sides=2

#for layerOffset in {0..62..1}
#do
    for PEthreshold in {6..40..2}
    do
        for PEYield in {10..68..2}  #PE per SiPM @ 1m away from SiPM of 3m long counter
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
#done
