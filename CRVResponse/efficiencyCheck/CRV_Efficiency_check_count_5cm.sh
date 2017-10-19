#layerOffset=42
timeWindow=10
sides=2
PEthreshold=18
moduleGap=5

    for photonYield in {3189,3781,4359}  # 24,28,32 PE/SiPM 1m away from SiPM for 5cm wide / 3m long counter
    do
      for layerOffset in {0..62..1}
      do

              directory=/pnfs/mu2e/scratch/users/ehrlich/workflow
#              name=CRV_Efficiency_check_5cm
#              name=CRV_Efficiency_check_5cm10
#              name=CRV_Efficiency_check_5cm0_5000r
#              name=CRV_Efficiency_check_5cm0_6000r
#              name=CRV_Efficiency_check_5cm0_6000
#              name=CRV_Efficiency_check_5cm0_upstreamPlanes
#              name=CRV_Efficiency_check_5cm_5000r
#              name=CRV_Efficiency_check_5cm_6000r
#              name=CRV_Efficiency_check_5cm_6000
#              name=CRV_Efficiency_check_5cm_upstreamPlanes
#              name=CRV_Efficiency_check_5cm_6600
#              name=CRV_Efficiency_check_5cm10_6600
#              name=CRVa_Efficiency_check_5cm_verticalPlanes
              name=CRV_Efficiency_check_5cm_downstreamPlanes
#              name=CRV_Efficiency_check_5cm10_upstreamPlanes

              files=$directory/$name/outstage/'*.'$name'_'gap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield/*/*/log.*.log
#              files=`ls $directory/$name/outstage/'*.'$name'_'gap$moduleGap'_'layerOffset$layerOffset'_'photonYield$photonYield/*/*/log.*.log`

              events=0
              eventsCoincidence=0

              for file in $files;
              do

                searchString="SUMMARY CrvCoincidencePE$PEthreshold""T$timeWindow"
#                searchCommand="grep '$searchString' $file"
                searchCommand="tail -n 150 $file | grep '$searchString'"
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
