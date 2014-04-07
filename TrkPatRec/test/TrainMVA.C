// This script trains the MVAs used to flag background (Compton)
// electron hits in the TTracker.  It processes the root files created by
// running the track reconstruction sequence with the following fcl parameter set:
// physics.producers.FlagBkgHits.diagLevel : 2
// We recommend processing a minimum of 10K events, fully mixed with backgrounds,
// as training sample.
// To run this script from the command line, invoke as:
// root -b -q "TrkPatRec/test/TrainMVA.C(\"files.txt\",2)" 
// where 'files.txt' should include the names of all the TFiles created by
// All 4 MVAs need to be trained.  By default, output files will go to the 
// 'tmva' directory in your offline release.
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  !!!!! MAKE SURE THAT DIRECTORY EXISTS BEFORE INVOKING THIS SCRIPT!!!!
// !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//  Those files can be copied into the
// correct location using the script TrkPatRec/test/InstallMVAWeights.sh
// Created by David Ding and David Brown, LBNL (April, 2014)

#include <iostream>

enum training{none=0,stereohit,nonstereohit,stereocluster,nonstereocluster};
void TrainMVA(const char* files="deltafiles.txt",training itrain=stereohit) {
  // fill TChain from files
  TChain* mytree = new TChain("FlagBkgHits/ddiag");
  std::ifstream fs(files,std::ifstream::in);
  char file[100];
  fs.getline(file,100);
  while(fs.good()){
    std::cout << "Adding file " << file << " to training sample" << std::endl;
    mytree->Add(file);
    fs.getline(file,100);
  }
  if(itrain == none){
    std::cout << "Test run, no training executed" << std::endl;
  } else if(itrain == stereohit) {
    std::cout << "Training Stereo Hits" << std::endl;
    gROOT->LoadMacro("TrkPatRec/test/TrainStereoHitsMVA.C+");
    TrainStereoHitsMVA(mytree);
  } else if (itrain == nonstereohit) {
    std::cout << "Training NonStereo Hits" << std::endl;
    gROOT->LoadMacro("TrkPatRec/test/TrainNonStereoHitsMVA.C+");
    TrainNonStereoHitsMVA(mytree);
  } else if(itrain == stereocluster) {
    std::cout << "Training Stereo Clusters" << std::endl;
    gROOT->LoadMacro("TrkPatRec/test/TrainStereoClusterMVA.C+");
    TrainStereoClusterMVA(mytree);
  } else if (itrain == nonstereocluster) {
    std::cout << "Training NonStereo Clusters" << std::endl;
    gROOT->LoadMacro("TrkPatRec/test/TrainNonStereoClusterMVA.C+");
    TrainNonStereoClusterMVA(mytree);
  }
  exit();
}

