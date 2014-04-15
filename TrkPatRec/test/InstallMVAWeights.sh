#!/bin/bash
baseDir="tmva"
while getopts "o:" opt ; do 
  case $opt in
    o)
      baseDir=$OPTARG
    ;;
    \?)
      echo "Invalid option, aborting"
      exit 1
    ;;
  esac
done
echo "Copying files from directory $baseDir"
cp $baseDir/TrainStereoHits/TMVAClassification_MLP.weights.xml TrkPatRec/test/StereoHits.weights.xml
cp $baseDir/TrainNonStereoHits/TMVAClassification_MLP.weights.xml TrkPatRec/test/NonStereoHits.weights.xml
cp $baseDir/TrainStereoCluster/TMVAClassification_MLP.weights.xml TrkPatRec/test/StereoCluster.weights.xml
cp $baseDir/TrainNonStereoCluster/TMVAClassification_MLP.weights.xml TrkPatRec/test/NonStereoCluster.weights.xml

