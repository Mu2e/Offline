#To submit the jobs:
#setup gridexport
#Delete content of CRVResponse/efficiency/submit
#Run this script from workspace/Offline
#Run gridexport from workspace/Offline and copy the location of the zip file to the joblist script to $CODELOCATION
#Run the joblist script from workspace/Offline

i=0
name=CRV_Efficiency_check_5cm_6000
layerOffset=42
moduleGap="3.0"
gapSmall="0.1"
gapLarge="0.5"
joblist=CRVResponse/efficiencyCheck/jobs.sh
#  for moduleGap in {"4.5","6.5"}
#  do
#    for layerOffset in {0..62..2}
#    do
      for PEYield in {32..68..2}  #PE per SiPM @ 1m away from SiPM of 3m long counter
      do

        ((i++));

        photonYield=`echo "47000.0*$PEYield/68.0" | bc -l`
        dz=`echo "800.0+$moduleGap+7*$gapLarge+8*$gapSmall" | bc -l`

        genconfigfile=CRVResponse/efficiencyCheck/submit/genconfig'_'$name'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/genconfig_5cm_6000.txt\"" >| $genconfigfile
        echo "double cosmicFromTH2.dz = $dz;" >> $genconfigfile

        geomfile=CRVResponse/efficiencyCheck/submit/geom'_'$name'_'$i.txt
        echo "#include \"CRVResponse/efficiencyCheck/geom_5cm_6000.txt\"" >| $geomfile
        echo "double crs.gapBetweenModules = $moduleGap;" >> $geomfile
        echo "double crs.gapSmall          = $gapSmall;" >> $geomfile
        echo "double crs.gapLarge          = $gapLarge;" >> $geomfile
        echo "double crs.layerOffset       = $layerOffset;" >> $geomfile
#        echo "vector<double> crs.offsetDirectionTest = {0, 0, 1};" >> $geomfile   #reverse layer offset

        fclfile=CRVResponse/efficiencyCheck/submit/$name'_'$i.fcl
        echo "#include \"CRVResponse/efficiencyCheck/CRV_Efficiency_check_5cm_6000.fcl\"" >| $fclfile
        echo "services.GeometryService.inputFile                      : \"$geomfile\"" >> $fclfile
        echo "physics.producers.generate.inputfile                    : \"$genconfigfile\"" >> $fclfile
        echo "physics.producers.CrvPhotons.scintillationYield         : $photonYield" >> $fclfile

        generate_fcl --description=$name --dsconf=$i --run=1 --events=50000 --njobs=20 $fclfile
        tar -zcvf fcls.$name.$i.tar.gz 000/cnf.$USER.$name.$i.*.fcl
        rm -rf 000
        mv fcls.$name.$i.tar.gz /pnfs/mu2e/scratch/outstage/ehrlich/fcls/.
        clustername=$name'_'moduleGap$moduleGap'_'gapLarge$gapLarge'_'gapSmall$gapSmall'_'layerOffset$layerOffset'_'PEYieldYield$PEYield
        echo "mu2eprodsys --expected-lifetime=4h --code=\$CODELOCATION --fcllist=/pnfs/mu2e/scratch/outstage/ehrlich/fcls/fcls.$name.$i.tar.gz --clustername=$clustername --dsconf=$i --wfproject=$name" >> $joblist

      done
#    done
#  done
