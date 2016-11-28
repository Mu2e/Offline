// This script trains the MVAs used to flag background (Compton)
// electron hits in the TTracker.  It processes the root files created by
// running the track reconstruction sequence with the following fcl parameter set:
// physics.producers.FlagBkgHits.diagLevel : 1
// physics.producers.FlagBkgHits.StereoHitMVACut : -1.0
// physics.producers.FlagBkgHits.NonStereoHitMVACut : -1.0
// then, for cluster training, rerun with:
// physics.producers.FlagBkgHits.StereoHitMVACut : 'good value'
// physics.producers.FlagBkgHits.NonStereoHitMVACut : 'good value'
// physics.producers.FlagBkgHits.StereoHitMVA.MVAWeights : 'tmva'.TrainStereoHits.TMVAClassification_MLP.weights.xml 
// physics.producers.FlagBkgHits.NonStereoHitMVA.MVAWeights : 'tmva'.TrainNonStereoHits.TMVAClassification_MLP.weights.xml 
//
// physics.producers.FlagBkgHits.StereoClusterMVACut : -1.0
// physics.producers.FlagBkgHits.NonStereoClusterMVACut : -1.0
// We recommend processing 10K events, fully mixed with backgrounds,
// as training sample.
// As an exanmple, to run this script from the command line to train stereo clusters, invoke as:
//  root -b -q "TrkPatRec/test/TrainMVA.C+(\"files.txt\",3,\"tmva\")"
// where 'files.txt' should include the names of all the TFiles created by FlagBkgHits
// when run on mixed Conversion Electron + background frames, with diagLevel >=2.
// All files will go to appropriately named subdirectories in the
// 'tmva' directory in your offline release.  The weight files can be copied into the
// correct location using the script TrkPatRec/test/InstallMVAWeights.sh
// Created by David Ding and David Brown, LBNL (April, 2014)

#include <iostream>
#include <fstream>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include "TrkPatRec/test/CreateMVAFactory.h"
#include "TChain.h"
#include "TCut.h"
#include "TFile.h"

enum training{none=0,stereohit,nonstereohit,stereocluster,nonstereocluster};
void TrainMVA(const char* files="deltafiles.txt",training itrain=stereohit,const char* basedir="tmva") {
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
//selection cuts
  TCut goodclu("tmean>700 && nprimary/nchits > 0.8 && nchits>4 && ns>1");
  TCut signalclu("pproc<20"); // signal are electrons from electromagnetic processes (Compton, delta, ...)
  TCut bkgclu("pgen==2&&mcmom>100"); // background are conversion electron hits
  TCut stclu("ngdstereo/ngdhits >0.49");
  TCut nonstclu("ngdstereo/ngdhits <0.49");
  TCut goodhit("pproc<20"); // only consider clusters from deltas for hit training
  TCut signalhit("_relation>=0"); // include all direct relations as signal
  TCut bkghit("_relation<0"); // unrelated hits to parent
  TCut sthit("_stereo");
  TCut nonsthit("!_stereo");
  TCut signalCut, backgrCut;
  
  std::vector<std::string> varnames;
  std::vector<std::string> vardescrip;
  std::string baseDir(basedir);
  std::string outDir;
  std::string weightexp;
  if(itrain == none){
    std::cout << "Test run, no training executed" << std::endl;
    exit(11);
  } else if(itrain == stereohit || itrain == nonstereohit) {
    varnames.push_back("_dphi");
    varnames.push_back("_drho");
    varnames.push_back("_dt");
    vardescrip.push_back("Azimuthal Angular Distance"); 
    vardescrip.push_back("Radial Distance");
    vardescrip.push_back("Temporal Distance");
    if(itrain == stereohit) {
      std::cout << "Training Stereo Hits" << std::endl;
      signalCut = goodclu+goodhit+sthit+signalhit;
      backgrCut = goodclu+goodhit+sthit+bkghit;
      outDir = baseDir + std::string("/TrainStereoHits");
    } else {
      std::cout << "Training NonStereo Hits" << std::endl;
      signalCut = goodclu+goodhit+nonsthit+signalhit;
      backgrCut = goodclu+goodhit+nonsthit+bkghit;
      outDir = baseDir + std::string("/TrainNonStereoHits");
    }
  } else if(itrain == stereocluster || itrain == nonstereocluster) {
    varnames.push_back("prho");
    varnames.push_back("srho");
    varnames.push_back("zmin");
    varnames.push_back("zmax");
    varnames.push_back("zgap");
    varnames.push_back("ns");
    varnames.push_back("nsmiss");
    varnames.push_back("sphi");
    varnames.push_back("ngdhits");
    vardescrip.push_back("Transverse radius");
    vardescrip.push_back("rho RMS"); 
    vardescrip.push_back("z of first hit");
    vardescrip.push_back("z of last hit");
    vardescrip.push_back("Biggest z gap between hits");
    vardescrip.push_back("Number of Stations");
    vardescrip.push_back("Number of Stations Missed");
    vardescrip.push_back("Phi RMS");
    vardescrip.push_back("Number of good hits");
    weightexp = "nprimary"; // weight each cluster by the # of primary hits
    if(itrain == stereocluster) {
      std::cout << "Training Stereo Clusters" << std::endl;
      signalCut = goodclu+stclu+signalclu;
      backgrCut = goodclu+stclu+bkgclu;
      outDir = baseDir + std::string("/TrainStereoCluster");
    } else { 
      std::cout << "Training NonStereo Clusters" << std::endl;
      signalCut = goodclu+nonstclu+signalclu;
      backgrCut = goodclu+nonstclu+bkgclu;
      outDir = baseDir + std::string("/TrainNonStereoCluster");
    }
  }
  // check base dir
  int stat = mkdir(baseDir.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (!stat)
    std::cout << "Directory " << baseDir << " created " << std::endl;
  stat = mkdir(outDir.c_str(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  if (!stat)
    std::cout << "Directory " << outDir << " created " << std::endl;
  // configure TMVA output: too bad this is a global!!!
  (TMVA::gConfig().GetIONames()).fWeightFileDir = outDir.c_str();
// create the factory
  std::string rootfile = outDir + std::string("/TMVAoutput.root");
  TFile* rootFile = TFile::Open( rootfile.c_str(), "RECREATE" );
  std::cout << "Defining signal with cut: " << signalCut.GetTitle() << std::endl;
  std::cout << "Defining background with cut: " << backgrCut.GetTitle() << std::endl;
  TMVA::Factory* factory = CreateMVAFactory(mytree,varnames, vardescrip,
  signalCut,backgrCut,rootFile,weightexp);
// run the training
  factory->TrainAllMethods();
  // ---- Evaluate all MVAs using the set of test events
  factory->TestAllMethods();
  // ----- Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();
  // Save the output
  rootFile->Close();
  std::cout << "==> Wrote root file: " << rootFile->GetName() << std::endl;
  std::cout << "==> TMVA Output in " << TMVA::gConfig().GetIONames().fWeightFileDir << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;
  delete factory;
// done!
  exit(0);
}

