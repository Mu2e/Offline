
// FMtool.cpp
 
#include "FigureOfMerit/inc/FMtool.h"

// C++ includes.
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <cassert>

// art and art externals includes
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/intermediate_table.h"
#include "fhiclcpp/make_ParameterSet.h"
#include "fhiclcpp/parse.h"

// mu2e includes
#include "FigureOfMerit/inc/PoissonCDFsolver.h"

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TTree.h"

namespace mu2e {



} // End of namespace mu2e for helper functions


using namespace std;
using fhicl::ParameterSet;

namespace mu2e {

FMtool::FMtool(ParameterSet const & p, std::ostream & o) 
  : pset( p )
  , os  ( o ) 
  , RPCtimeProbability(RPCtimeProbabilitySpline())
{
  decideVerbosity();  // under control of the fcl file
}

void FMtool::analyze()
{
  // TODO - replace the setup steps by a driver that uses a file or whatever
  //        to get data inputs.
  
  setFilterLevels();
        if (OUTPUT_detailedTrace) os << "### setFilterLevels returned \n";
  setCanonicals();
        if (OUTPUT_detailedTrace)  os << "### setCanonicals returned \n";
  setRootGraphics();
        if (OUTPUT_detailedTrace) os << "### setRootGraphics returned \n";  
	if (rg.enabled) {
	  establishHistograms();
	  if (OUTPUT_detailedTrace) os << "### establishHistograms returned \n";
	}  
  double nCE = obtainCEdata();
        if (OUTPUT_detailedTrace)  os << "### obtainCEdata returned \n";
  normalizeSignalEfficiency(nCE);  
        if (OUTPUT_detailedTrace)  os << "### normalizeSignalEfficiency returned \n";
  bool use_diowt;
  double genDIO = obtainDIOdata(use_diowt);
        if (OUTPUT_detailedTrace)  os << "### obtainDIOdata returned \n";
  normalizeDIObackground(genDIO,use_diowt);  
        if (OUTPUT_detailedTrace)  os << "### normalizeDIObackground returned \n";
  double genRPC = obtainRPCdata();
        if (OUTPUT_detailedTrace)  os << "### obtainRPCdata returned \n";
  normalizeRPCbackground(genRPC);  
        if (OUTPUT_detailedTrace)  os << "### normalizeRPCbackground returned \n";
  applyFofM();
        if (OUTPUT_detailedTrace)  os << "### applyFofM returned \n";
  closeRootFiles();
        if (OUTPUT_detailedTrace)  os << "### closeRootFiles returned \n";
}

void FMtool::setFilterLevels()
{
  ParameterSet cuts = pset.get<ParameterSet>("trackerCuts");

  minimum_nactive   = cuts.get<int>  ("minimum_nactive"); 
  maximum_t0err     = cuts.get<float>("maximum_t0err"); 
  maximum_fitmomerr = cuts.get<float>("maximum_fitmomerr"); 
  minimum_fitcon    = cuts.get<float>("minimum_fitcon"); 
  tCuts = pset.get< std::vector<double> >("t0cuts");

  recordFclGroup("trackerCuts");
  recordFcl("minimum_nactive  ",minimum_nactive);
  recordFcl("maximum_t0err    ",maximum_t0err);
  recordFcl("maximum_fitmomerr",maximum_fitmomerr);
  recordFcl("minimum_fitcon   ",minimum_fitcon);

  for (size_t i = 0; i < tCuts.size(); ++i) {
    recordFcl("tCuts          ",tCuts[i]);
  }
}  // setFilterLevels


  void FMtool::setCanonicals() 
{
  // TODO - allow this to be set from file or whatever.

  // Canonical range
  canonicalRangeLo = pset.get<double>("canonicalRangeLo",103.5);
  canonicalRangeHi = pset.get<double>("canonicalRangeHi",104.7);
  recordFclGroup("Canonical range for early counts of inputs");
  recordFcl("canonicalRangeLo", canonicalRangeLo);
  recordFcl("canonicalRangeHi", canonicalRangeHi);
  
  ParameterSet binningPset = pset.get<ParameterSet>("binning",ParameterSet());
  double lowestBinDefault = 102.5;  
  double topOfLastBinDefault = 106.0;  
  stoppedMuonsDef = 143107;
  stoppedMuonsThisRun = 143107;
  sc_factor = 1.0;
  int nBinsDefault = 350;
  lowestBin    = binningPset.get<double>("lowestBin",   lowestBinDefault);
  topOfLastBin = binningPset.get<double>("topOfLastBin",topOfLastBinDefault);
  nBins        = binningPset.get<int>   ("nBins",       nBinsDefault);
  nSplineBins  = binningPset.get<int>   ("nSplineBins", nBins);
  stoppedMuonsDef  = binningPset.get<int>   ("stoppedMuonsDef", stoppedMuonsDef);
  stoppedMuonsThisRun  = binningPset.get<int>   ("stoppedMuonsThisRun", stoppedMuonsThisRun);

  binSize = (topOfLastBin - lowestBin)/nBins;
  recordFclGroup("binning");
  recordFcl("lowestBin   ", lowestBin);
  recordFcl("topOfLastBin", topOfLastBin);
  recordFcl("nBins       ", nBins);
  recordFcl("nSplineBins ", nSplineBins);
  recordFcl("stoppedMuonsDef ", stoppedMuonsDef);
  recordFcl("stoppedMuonsThisRun ", stoppedMuonsThisRun);  

  // Experiment quantities
    // TODO -- obtain these normalization-related constants from some sort of 
    // database or conditions service. 
  protonsOnTarget = 3.6e20;  
  cadence = 1694;  // nsec after start of proton pulse
  sc_factor = double(stoppedMuonsThisRun)/double(stoppedMuonsDef);  //this facor is obtained w.r.t. def value
  stoppedMuonsPerPOTold = 0.0016;  // I have seen .00215; new  number is .0016
  stoppedMuonsPerPOT = stoppedMuonsPerPOTold*sc_factor;
  //std::cout<<"     scale factor: "<<sc_factor<<"   "<<stoppedMuonsPerPOT<<endl;
  capturedMuonsPerStoppedMuon = 0.609; // DocDB 48 - I have seen .59
  RPCperStoppedPion = 0.021;     // DocDB 1087
  stoppedPionsPerPOT = 1.53e-6;

  recordFclGroup("physics (hardwired)");
  recordFcl("protonsOnTarget   ", protonsOnTarget);
  recordFcl("cadence           ", cadence);
  recordFcl("RPCperStoppedPion ", RPCperStoppedPion);
  recordFcl("stoppedPionsPerPOT ", stoppedPionsPerPOT);
  recordFcl("stoppedMuonsPerPOT without scale factor",  stoppedMuonsPerPOTold);
  recordFcl("stoppedMuons default case (expected num 143107 for 17 foils)",  stoppedMuonsDef);
  recordFcl("stoppedMuons for this job",  stoppedMuonsThisRun);
  recordFcl("scale factor (stoppedMuonsDef/ThisRun) multipied to stoppedMuonsPerPOT", sc_factor);
  recordFcl("used stoppedMuonsPerPOT for default case: ", stoppedMuonsPerPOTold);
  recordFcl("used stoppedMuonsPerPOT for this run: ", stoppedMuonsPerPOT);

  // control of table size and returned vectors from FofM
  maximumSignalCount = 25;
 
} // setCanonicals
 
void FMtool::setRootGraphics()
{
  // int branchStyle = 0;  // old branch style

  rg.rootFile = 
        pset.get<std::string>("Root.rootFile","noRootFileNameSpecified");
  if (rg.rootFile != "noRootFileNameSpecified") {
    rg.enabled = true;
    std::string cintFileDefault = rg.rootFile + ".cint"; 
    rg.cintFile = pset.get<std::string>("Root.cintFile",cintFileDefault);
    rg.rootfp = new TFile (rg.rootFile.c_str(),"RECREATE","FMtool generated ROOT file");
    rg.tree = new TTree("treeout","for stored values");
    int ivalv = 0;
    rg.tree->Branch("timeCut", &rg.valv[ivalv], "timeCut/F");    ++ivalv; 
    rg.tree->Branch("momWLow", &rg.valv[ivalv], "momWLow/F");    ++ivalv; 
    rg.tree->Branch("momWHig", &rg.valv[ivalv], "momWHig/F");    ++ivalv; 
    rg.tree->Branch("90PCL",   &rg.valv[ivalv], "90PCL/F");      ++ivalv; 
    rg.tree->Branch("minMeaningfullDIOtail",   &rg.valv[ivalv], "minMeaningfullDIOtail/F");      ++ivalv; 
    rg.tree->Branch("canonicalRangeLo",   &rg.valv[ivalv], "canonicalRangeLo/F");      ++ivalv; 
    rg.tree->Branch("canonicalRangeHi",   &rg.valv[ivalv], "canonicalRangeHi/F");      ++ivalv; 
    rg.tree->Branch("lowestBin",   &rg.valv[ivalv], "lowestBin/F");      ++ivalv; 
    rg.tree->Branch("topOfLastBin",   &rg.valv[ivalv], "topOfLastBin/F");      ++ivalv; 
    rg.tree->Branch("stoppedMuonsDef",   &rg.valv[ivalv], "stoppedMuonsDef/F");      ++ivalv; 
    rg.tree->Branch("stoppedMuonsThisRun",   &rg.valv[ivalv], "stoppedMuonsThisRun/F");      ++ivalv; 
    rg.tree->Branch("sc_factor",   &rg.valv[ivalv], "sc_factor/F");      ++ivalv; 
    rg.tree->Branch("stoppedMuonsPerPOTold",   &rg.valv[ivalv], "stoppedMuonsPerPOTold/F");      ++ivalv; 
    rg.tree->Branch("stoppedMuonsPerPOT",   &rg.valv[ivalv], "stoppedMuonsPerPOT/F");      ++ivalv; 
                                                                 
    if (rg.rootfp->IsZombie()) {
      os << "Failed to properly open the root output file " 
                << rg.rootFile << "\n";
      std::cerr << "Failed to properly open the root output file " 
                << rg.rootFile << "\n";
      delete rg.rootfp;
      std::exit(1);
    }
    os << "Root file to be created: " << rg.rootFile
       << "\ncint file to be created: " << rg.cintFile << "\n";
    // TTree::SetBranchStyle(branchStyle);
    // rg.tree = new TTree ("input","input data");
    // rg.rootfp->Add(rg.tree);  // The examples do not say this is needed, but...
    // rg.tree->SetDirectory(rg.rootfp); // The examples do not say this is needed,
                                      // RootOutputTree.cc does it...
  } else {
    rg.enabled = false;
    os << "No fcl parameter Root.rootFile detected. \n"
       << "No root file containing data for graphics will be produced.\n";
    std::cerr << "No fcl parameter Root.rootFile detected. \n"
       << "No root file containing data for graphics will be produced.\n";
  }

}  // setRootGraphics

void FMtool::establishHistograms() 
{
  if (rg.enabled) {

    Long_t nn=0;   const int nhist = rg.nhist;
    double xmin[nhist], xmax[nhist], nsiz[nhist];  int nbins[nhist]; TString htil;    TString Hnn;  

    for(int i=0; i<nhist; ++i){
      rg.h_fitmomCE[i]=NULL;          rg.h_fitmomDIO[i]=NULL;
      rg.hQ_fitmomCE[i]=NULL;         rg.hQ_fitmomDIO[i]=NULL;
    }
    
    int ibin=0;
    nsiz[ibin]=0.5;   xmin[ibin]=-1.0;   xmax[ibin]=149.0;    ++ibin;
    nsiz[ibin]=0.1;   xmin[ibin]=96.0;   xmax[ibin]=106.0;    ++ibin;      
    nsiz[ibin]=0.05;   xmin[ibin]=103.0;   xmax[ibin]=106.0;  ++ibin;

    if(ibin!=nhist)cout<<"Error in defining hist "<<ibin<<"   "<<nhist<<endl;
    
    for(int i=0; i<nhist; ++i)nbins[i] = int ( TMath::Abs(xmax[i]-xmin[i])/nsiz[i]);
    
    for(int i=0; i<nhist; ++i){
      nn=i;  
      //cout<<i<<" --> Hist    xmin/max/nbins "<<xmin[i]<<"  "<<xmax[i]<<"   "<<nbins[i]<<endl;
      Hnn="h_fitmomCE"; htil="CE fitmom";  Hnn=Hnn+nn;  rg.h_fitmomCE[i]=new TH1F(Hnn, htil, nbins[i], xmin[i], xmax[i]);
      Hnn="h_fitmomDIO"; htil="DIO fitmom";  Hnn=Hnn+nn;  rg.h_fitmomDIO[i]=new TH1F(Hnn, htil, nbins[i], xmin[i], xmax[i]);
      Hnn="hQ_fitmomCE"; htil="CE fitmom, Trk Qual";  Hnn=Hnn+nn;  rg.hQ_fitmomCE[i]=new TH1F(Hnn, htil, nbins[i], xmin[i], xmax[i]);
      Hnn="hQ_fitmomDIO"; htil="DIO fitmom, Trk Qual";  Hnn=Hnn+nn;  rg.hQ_fitmomDIO[i]=new TH1F(Hnn, htil, nbins[i], xmin[i], xmax[i]);
   }

    rg.hCEspectrum = new TH1F("CEfitmomX", "CE fitted momentum", 100, 96.0, 106.0);
    rg.hCEspectrumX = new TH1F("CEfitmom", "CE fitted momentum", nBins, lowestBin, topOfLastBin);

    rg.hDIOspectrum = new TH1F("DIOfitmomX", "DIO fitted momentum", 100, 96.0, 106.0);
    rg.hDIOspectrumX = new TH1F("DIOfitmom", "DIO fitted momentum", nBins, lowestBin, topOfLastBin);

     //rg.cint <<  "// Established histogram DIOfitmom -  DIO tracks passing cuts \n"
     //        <<  "//                                    within the binning area \n";
     //rg.hDIOspectrumX = new TH1F("DIOfitmomX", "DIO fitted momentum", 100, 96.0, 106.0);
     //rg.cint <<  "// Established histogram DIOfitmomX - DIO tracks passing cuts \n"
     //        <<  "//                                    with p from 96 to 106 MeV \n";
  }
} // establishHistograms

bool FMtool::RootGraphics::write_cint_file() const 
{
  if (enabled) {
    std::ofstream ofs (cintFile.c_str());
    if (!ofs) {
      std::cerr << "Failed to open cint file " << cintFile << " for writing\n";
      return false;
    }
    ofs << cint.str();
  }
  return true;
}

void FMtool::closeRootFiles() 
{
  if (rg.enabled) {  

    //rg.valsD.Write("valsD");  

    //std::cerr << "### closeRootFiles entered \n"; 
    if (OUTPUT_detailedTrace)  os << "### closeRootFiles entered \n";
    rg.write_cint_file();
    //std::cerr << "### write_cint_file returned \n"; 
    if (OUTPUT_detailedTrace)  os << "### write_cint_file returned \n";
    rg.rootfp->cd();   
    rg.rootfp->Write();                                 
    rg.rootfp->Close();  
  }
  //std::cerr << "### closeRootFiles returned \n"; 
  if (OUTPUT_detailedTrace)  os << "### closeRootFiles returned \n";
}

double FMtool::obtainCEdata() 
{
  // Obtain the CE data as a vector of momenta

  // Returns a number of CE's generated, for use when normalizing.

  if (OUTPUT_progress || OUTPUT_signalEfficiency) {
    os << "\nObtaining conversion electron signal efficiency spectrum...\n";
  }
  
  ParameterSet CEpset = pset.get<ParameterSet>("conversionElectrons");
  std::vector<std::string> CEfileList =
    CEpset.get< std::vector<std::string> > ("CEdataFiles");

  double numberOfGeneratedMuonConversions = 
    CEpset.get<double>("numberOfGeneratedMuonConversions");
 
  CEliveGateFraction = CEpset.get<double>("liveGateFraction", 1.0);
  adHocSignalrescaler = CEpset.get<double>("adHocSignalrescaler",1.0); 
        recordFclGroup("conversionElectrons");
        recordFcl("numberOfGeneratedMuonConversions", numberOfGeneratedMuonConversions); 
        recordFcl("CEliveGateFraction              ",  CEliveGateFraction);
        if (adHocSignalrescaler != 1.0) {
          recordFcl("adHocSignalrescaler             ",  adHocSignalrescaler);
        }
        recordFclGroup("conversionElectrons.CEdataFiles");
	for (size_t i=0; i< CEfileList.size(); ++i) {
	  recordFcl("", CEfileList[i]);
	}
  size_t nTcuts = tCuts.size();
  sigEfficiency.clear();
  sigEfficiency.resize(nTcuts);
  for (size_t tCutNumber = 0; tCutNumber < nTcuts; ++tCutNumber) {
    sigEfficiency[tCutNumber].clear();
    sigEfficiency[tCutNumber].resize(nBins);
  }
  extractFitmom ( CEfileList, sigEfficiency );

  extractFitmom ( CEfileList, sigEfficiency, false, false, true ); 

  if (OUTPUT_dataProperties) {
    os << "CE's based on " << numberOfGeneratedMuonConversions 
       << " generated muon conversions \n"; 
    if (CEliveGateFraction != 1.0) {
      os << "CE data used a live time gate fraction of " 
         << CEliveGateFraction << "\n";
    }
  }

  return numberOfGeneratedMuonConversions;
  
} // obtainCEdata

void FMtool::normalizeSignalEfficiency(double numberOfGeneratedMuonConversions)
{
  // Normalization of the signal efficiencies
  // Physics Notes: 
  //    *  Our convention is to normalize the signal efficiency to a
  //       per POT number.
  //    *  Only muons that are captured by the nucleus are 
  //       candidates for direct conversion.  Thus we need to divide
  //       by capturedMuonsPerStoppedMuon to go from number of CE's generated
  //       to number of POT's represented.

  liveGateFraction = CEliveGateFraction;

  double POTsRepresented = numberOfGeneratedMuonConversions /
     ( stoppedMuonsPerPOT * capturedMuonsPerStoppedMuon * liveGateFraction );
  double sigEffCountNormalization =  1.0/POTsRepresented;
  double sigEffNormalization =  sigEffCountNormalization / binSize;

  os << "|| stoppedMuonsPerPOT          s = " << stoppedMuonsPerPOT << "\n";
  os << "|| capturedMuonsPerStoppedMuon c = " << capturedMuonsPerStoppedMuon << "\n";
  os << "|| CEliveGateFraction          F = " << liveGateFraction << "\n";
  os << "||                         s*c*F = " << stoppedMuonsPerPOT*
                          capturedMuonsPerStoppedMuon*liveGateFraction << "\n";
  os << "|| Divided by numberOfGeneratedMuonConversions (" << numberOfGeneratedMuonConversions
     << ") = " <<   stoppedMuonsPerPOT*
                          capturedMuonsPerStoppedMuon*liveGateFraction 
                          / numberOfGeneratedMuonConversions << "\n"; 
  os << "|| sigEffCountNormalization = " 
     << sigEffCountNormalization << "\n";                   
  os << "|| Times Luminosity " << protonsOnTarget << " = " 
     << sigEffCountNormalization * protonsOnTarget << "\n";             
  os << "|| (This is number of CE counts in experiment per CE in the input files, \n"
     << "||  *if* the branching ratio were 1) \n";

  if (OUTPUT_signalEfficiency) {
    os << "sigEffNormalization (count normalization/bin size) = " 
       << sigEffNormalization << "\n";
  }


  if ( adHocSignalrescaler != 1.0 ) {
    os << "\n** For this trial, signal normalization is scaled to "
              << adHocSignalrescaler << " times its actual value **\n\n";
    sigEffNormalization *= adHocSignalrescaler;
  }

  for (size_t t = 0; t < tCuts.size(); ++t) {
    if (OUTPUT_signalEfficiency)
      { os << "Signal Efficiency with time cut at " << tCuts[t] << ": \n"; }
    unsigned int canonicalRangeCEcount = 0;
    for (unsigned int i=0; i < sigEfficiency[t].size(); ++i) {
      double p = lowestBin + binSize*i;
      if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) ) {
        canonicalRangeCEcount += sigEfficiency[t][i];
      }
    }
    if (OUTPUT_signalEfficiency) {
      os << "There are " <<  canonicalRangeCEcount
         << " CE's in the files in the range of "
         << canonicalRangeLo << " - "  << canonicalRangeHi << " \n";
      os << "With normalization, this amounts to "
         << canonicalRangeCEcount * sigEffCountNormalization
                                  * protonsOnTarget * 1.0e-16
         << " CE events in range if BR = 1.0E-16\n";
    }

//    double priorIntegratedSignalEfficiencyInRange = 0;
//
//    for (unsigned int i=0; i < sigEfficiency[t].size(); ++i) {
//      double p = lowestBin + binSize*i ;
//      if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) ) {
//        priorIntegratedSignalEfficiencyInRange += sigEfficiency[t][i];
//      } 
//    }
//    os << "For time window " << tCuts[t] << ", before normalization, \n"
//      << "summed signal efficency in range (" << canonicalRangeLo
//       << ", " << canonicalRangeHi << ") is " 
//       << priorIntegratedSignalEfficiencyInRange << "\n";

    for (unsigned int i=0; i < sigEfficiency[t].size(); ++i) {
      sigEfficiency[t][i] *= sigEffNormalization;     // MULIPLICATION STEP AAA
    }

    if (OUTPUT_signalEfficiency && OUTPUT_spectra) {
      for (unsigned int i=0; i < sigEfficiency[t].size(); ++i) {
        os << lowestBin + binSize*i << ":   " 
           << sigEfficiency[t][i]/sigEffNormalization;
        os << " --> " << sigEfficiency[t][i] 
           << " ==> " << sigEfficiency[t][i] * protonsOnTarget << "\n";
      }
    }
    if (OUTPUT_signalEfficiency) {
      double integratedSignalEfficiencyInRange = 0;
      for (unsigned int i=0; i < sigEfficiency[t].size(); ++i) {
        double p = lowestBin + binSize*i ;
        if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) ) {
          integratedSignalEfficiencyInRange += sigEfficiency[t][i];
        } 
      }
      os << "Summed signal efficency in range (" << canonicalRangeLo
         << ", " << canonicalRangeHi << ") is " 
         << integratedSignalEfficiencyInRange << "\n";
    }   
  } // end of loop over t
  
} // normalizeSignalEfficiency

double FMtool::obtainDIOdata(bool & use_diowt) 
{
  // TODO - get the instructions for obtaining this data 
  // from a file or from the user.

  if (OUTPUT_progress || OUTPUT_backgroundStrength) {
    os << "\nObtaining DIO spectrum...\n";
  }

  // Obtain the DIO data as a vector of momenta 

  ParameterSet DIOpset = pset.get<ParameterSet>("decayInOrbitElectrons");
  std::vector<std::string> DIOfileList =
    DIOpset.get< std::vector<std::string> > ("DIOdataFiles");
        recordFclGroup("decayInOrbitElectrons.DIOdataFiles");
	for (size_t i=0; i< DIOfileList.size(); ++i) {
	  recordFcl("", DIOfileList[i]);
	}
  double numberOfGeneratedDIOs = DIOpset.get<double>("numberOfGeneratedDIOs");
  DIOliveGateFraction = DIOpset.get<double>("liveGateFraction");
  use_diowt = DIOpset.get<bool>("dioTracksWeighted");
  if (use_diowt) {
    DIOflatGenerationWindowLo = 
        DIOpset.get<double>("DIOflatGenerationWindowLo ", 0);
    DIOflatGenerationWindowHi = 
        DIOpset.get<double>("DIOflatGenerationWindowHi ", 0);
    if (OUTPUT_detailedTrace) 
      os << "### DIOflatGenerationWindowLo = " << DIOflatGenerationWindowLo
         << "\nDIOflatGenerationWindowHi = " << DIOflatGenerationWindowHi
         << "\n";
    if ( (DIOflatGenerationWindowLo == 0) ||
         (DIOflatGenerationWindowHi == 0) ) {
      os << "Since dioTracksWeighted = true, "
         << "DIOflatGenerationWindowLo and Hi are needed \n"
         << "FMtool will exit because it cannot normalize the DIO background\n";
      std::exit(1);
    }         
    if ( DIOflatGenerationWindowLo >= DIOflatGenerationWindowHi ) {
      os << "DIOflatGenerationWindowLo = " << DIOflatGenerationWindowLo
         << "\nDIOflatGenerationWindowHi = " << DIOflatGenerationWindowHi
         << "\nThis is not a suitable range for having generated DIOs\n";
      std::exit(1);
    }
  } 
  //std::cerr << "### 2\n";
  DIObackground.clear();
  DIObackground.resize(nBins);
  double count103 = 0.0;
  double count104 = 0.0;
  double countRatioThreshold = 5.0;  
  countMcmom ( DIOfileList, count103, count104, use_diowt ); 
  //std::cerr << "### 3\n";

  os << "The ratio of DIO's in 102.5-103.5 to 103.5-104.5 in this set is "
     << count103/count104 << " -- \n";

  os << "That represents " << count103 << " and " << count104 
     << " counts respectively\n";
     
  bool diowt_discrepancy = false;
  if (count104 <= 0) {
    if ( use_diowt ) {
      std::cerr << "Warning: Discrepancy detected involving use of diowt! \n"; 
      os << "There are no DIO's near 104 in this set.\n"
         << "Apparently, the DIO tracks were generated \n"
         << "according to the DIO momentum distribution " 
         << "and should not not be weighted.\n"
         << "But dioTracksWeighted (in the fcl) has value true.\n";
      os        << "Warning: Discrepancy detected involving use of diowt! \n"; 
      diowt_discrepancy = true;
    }
  } else if (count103/count104 > countRatioThreshold){
    if ( use_diowt ) {
      std::cerr << "Warning: Discrepancy detected involving use of diowt! \n"; 
      os << "The ratio of DIO's near 103 to near 104 in this set is "
         << count103/count104 << " -- \n"
         << "Apparently, the DIO tracks were generated \n"
         << "according to the DIO momentum distribution " 
         << "and should not be weighted\n"
         << "But dioTracksWeighted (in the fcl) has value true.\n";
      os        << "Warning: Discrepancy detected involving use of diowt! \n"; 
      diowt_discrepancy = true;
    }
  } else {
    if ( !use_diowt ) {
      std::cerr << "Warning: Discrepancy detected involving use of diowt! \n"; 
      os << "The ratio of DIO's near 103 to near 104 in this set is "
         << count103/count104 << " -- \n"
         << "Apparently, the DIO tracks were generated \n"
         << "according to a roughly flat momentum distribution "
         << "and need to be weighted.\n" 
         << "But dioTracksWeighted (in the fcl) has value false.\n";
      os        << "Warning: Discrepancy detected involving use of diowt! \n"; 
      diowt_discrepancy = true;
    }
  }
  //std::cerr << "### 4\n";
  if (diowt_discrepancy) {
    os << "Discontinuing FMtool due to discrepancy in use of DIO weights\n";
    std::cerr << "Discontinuing FMtool due to discrepancy in use of DIO weights\n";
    std::exit(1);
  }
  if ( use_diowt ) {
    fractionOfDIOsRepresented = 1.0;
  } else {
    fractionOfDIOsRepresented = 
      DIOpset.get<double>("fractionOfDIOsRepresented",-99);
      if ( fractionOfDIOsRepresented > 0) {
	  recordFcl("fractionOfDIOsRepresented", fractionOfDIOsRepresented);
        os << "fractionOfDIOsRepresented is set by parameter to " 
           << fractionOfDIOsRepresented << "\n";
      }
    if (fractionOfDIOsRepresented < 0) {
      double DIOgenerationLowerLimit = 
      DIOpset.get<double>("DIOgenerationLowerLimit");
      CzarneckiDIOspectrumFunction cz(1.0);
      fractionOfDIOsRepresented = cz.integrated(DIOgenerationLowerLimit);
      os << "DIO Generation Lower Limit = " << DIOgenerationLowerLimit
         << "  fraction of DIOs represented = " 
         << fractionOfDIOsRepresented << "\n";
      os << "Fraction at " << canonicalRangeLo << " would be " 
          << cz.integrated(canonicalRangeLo) << " (" 
          << fractionOfDIOsRepresented/cz.integrated(canonicalRangeLo)
          << ")\n";
    }
  }
  //std::cerr << "### 5\n";
  minimumMeaningfulDIOtail = DIOpset.get<int>("minimumMeaningfulDIOtail",5); 

  bool alsoReadWithNoT0cut = DIOpset.get<bool>("alsoReadWithNoT0cut",false);
  if (alsoReadWithNoT0cut) {
    os << "DIO's ignoring t0 cuts:\n: "; 
    tracksExtracted = 0;
    tracksBinned = 0;
    size_t nTcuts = tCuts.size();
    DIObackground.clear();
    DIObackground.resize(nTcuts);
    for (size_t tCutNumber = 0; tCutNumber < nTcuts; ++tCutNumber) {
      DIObackground[tCutNumber].clear();
      DIObackground[tCutNumber].resize(nBins);
    }
    extractFitmom ( DIOfileList, DIObackground, true, use_diowt, false );  
    os << "Total DIO tracks extracted passing cuts (except t0): " 
       << tracksExtracted << "\nof which " << tracksBinned << " were binned \n";
    os << "Now re-reading DIO's not ignoring t0 cuts:\n: "; 
  }
  
  tracksExtracted = 0;
  tracksBinned = 0;
  size_t nTcuts = tCuts.size();
  DIObackground.clear();
  DIObackground.resize(nTcuts);
  for (size_t tCutNumber = 0; tCutNumber < nTcuts; ++tCutNumber) {
    DIObackground[tCutNumber].clear();
    DIObackground[tCutNumber].resize(nBins);
  }
  extractFitmom ( DIOfileList, DIObackground, true, use_diowt, true );  
  os << "Total DIO tracks extracted passing cuts: " 
     << tracksExtracted << "\nof which " << tracksBinned << " were binned \n";
   
  if (OUTPUT_dataProperties) {
    os << "DIO's based on " << numberOfGeneratedDIOs 
       << "  generated DIO electrons \n"; 
    if (DIOliveGateFraction != 1.0) {
      os << "DIO data used a live time gate fraction of " 
         << DIOliveGateFraction << "\n";
    }
  }

  adHocDIOrescaler = DIOpset.get<double>("adHocDIOrescaler",1.0);

        recordFclGroup("decayInOrbitElectrons");
        recordFcl("numberOfGeneratedDIOs     ", numberOfGeneratedDIOs); 
        recordFcl("liveGateFraction          ", liveGateFraction); 
        recordFcl("dioTracksWeighted         ", use_diowt ); 
        if (use_diowt) {
          recordFcl("DIOflatGenerationWindowLo ", DIOflatGenerationWindowLo);
          recordFcl("DIOflatGenerationWindowHi ", DIOflatGenerationWindowHi);
	}
        recordFcl("minimumMeaningfulDIOtail  ", minimumMeaningfulDIOtail);
        if ( adHocDIOrescaler != 1.0 ) {
          recordFcl("adHocDIOrescaler          ", adHocDIOrescaler);
        }

  makeDIOreferenceTF1(tracksExtracted, DIOflatGenerationWindowLo);

  return numberOfGeneratedDIOs;
  
} // obtainDIOdata

double FMtool::obtainRPCdata() 
{
  // Obtain the RPC data as a vector of momenta

  if (OUTPUT_progress || OUTPUT_backgroundStrength) {
    os << "\nObtaining RPC spectrum...\n";
  }

  ParameterSet RPCpset = 
    pset.get<ParameterSet>("radiativePionCaptureElectrons");
  std::vector<std::string> RPCfileList =
    RPCpset.get< std::vector<std::string> > ("RPCdataFiles");
        recordFclGroup("radiativePionCaptureElectrons.RPCdataFiles");
	for (size_t i=0; i< RPCfileList.size(); ++i) {
	  recordFcl("", RPCfileList[i]);
	}
  double numberOfGeneratedRPCs = 
    RPCpset.get<double>("numberOfGeneratedRPCs");

  // There is never a meaningful minimum_t0 to cut on in the RPC data 
  // because the distribution drops so rapidly one would never create an honest 
  // one.  Instead, analytic weighting is always used.  So we use an option
  // of extractFitmom that turns off time cutting, and produce only one
  // entry in the collection of time-cut-dependent "histograms."

  size_t nTcuts = tCuts.size();
  RPCbackground.clear();
  RPCbackground.resize(nTcuts);
  for (size_t tCutNumber = 0; tCutNumber < nTcuts; ++tCutNumber) {
    RPCbackground[tCutNumber].clear();
    RPCbackground[tCutNumber].resize(nBins);
  }
  RPCbackground[0].clear();
  RPCbackground[0].resize(nBins);

  extractFitmom    ( RPCfileList, 
                     RPCbackground,
                     false, //  dio -- this is not for dio
                     false, //  use_diowt
                     false  //  apply_timeCuts
                   );

  extinction = RPCpset.get<double>("extinction");
  
  if (OUTPUT_dataProperties) {
    os << "RPC's based on " << numberOfGeneratedRPCs 
       << " generated RPC electrons from stopped pions \n"; 
    os << "Using extinction of " << extinction << "\n";
  }

  adHocRPCrescaler = RPCpset.get<double>("adHocRPCrescaler",1.0);
        
        recordFclGroup("radiativePionCaptureElectrons");
	recordFcl("numberOfGeneratedRPCs",numberOfGeneratedRPCs);
	recordFcl("extinction           ",extinction);
	if (adHocRPCrescaler != 1.0) {
          recordFcl("adHocRPCrescaler     ",adHocRPCrescaler);
	}  
  
  return numberOfGeneratedRPCs;
  
} // obtainRPCdata

void FMtool::extractFitmom 
  ( TTree * tracks, 
    std::vector< std::vector<double> > & counts, 
    std::vector<NthHighest> & highPs,
    bool dio,
    bool use_diowt, 
    bool apply_timeCuts )
{
  // Set up place to obtain the momentum and other data for immediate use 
  tracks->SetBranchStyle(0);
  
  tracks->SetBranchStatus("*",0);
  tracks->SetBranchStatus("nactive",1);
  tracks->SetBranchStatus("t0err",1);
  tracks->SetBranchStatus("fitmomerr",1);
  tracks->SetBranchStatus("fitcon",1);
  tracks->SetBranchStatus("fitstatus",1);
  tracks->SetBranchStatus("fitmom",1);
  tracks->SetBranchStatus("diowt",1);
  tracks->SetBranchStatus("t0",1);
  int   nactive = 0;   tracks->SetBranchAddress("nactive"  , &nactive);
  float t0err = 0;     tracks->SetBranchAddress("t0err"    , &t0err);    
  float fitmomerr = 0; tracks->SetBranchAddress("fitmomerr", &fitmomerr);    
  float fitcon = 0;    tracks->SetBranchAddress("fitcon"   , &fitcon);    
  int   fitstatus = 0; tracks->SetBranchAddress("fitstatus", &fitstatus);    
  float fitmom = 0;    tracks->SetBranchAddress("fitmom"   , &fitmom);    
  float diowt = 0;     tracks->SetBranchAddress("diowt"    ,  &diowt);    
  float t0 = 0;        tracks->SetBranchAddress("t0"       ,  &t0);    
  int ntracks = 0;
  //std::cerr << "### 6\n";
  double binSize = (topOfLastBin - lowestBin)/nBins;
  std::vector<int> nBinnedTracks(tCuts.size());
  std::vector<double> dioTotalWeight(tCuts.size());
  std::vector<double> dioTotalWeightInCanonicalRange (tCuts.size());
  size_t nt = tCuts.size();
  int treeEntries = tracks->GetEntries();
  if (OUTPUT_detailedTrace) os << "### treeEntries = " << treeEntries << "\n";
  int nfitstatus=0;
  int nnactive =0;
  int nfitcon=0;
  int nfitmomerr = 0;
  int nt0err = 0;
  int nfitmom=0;
  //std::cerr << "### 7\n";
  for (int i = 0; i < treeEntries; ++i) {
    int bytesRead = tracks->GetEntry(i);
    if (bytesRead == 0) break;
    ++ntracks;
    if (OUTPUT_detailedCutTrace) {
      accumulateCutStatistics(  fitstatus,  nactive,  fitcon,  fitmomerr,  t0err, fitmom, 
			       nfitstatus, nnactive, nfitcon, nfitmomerr, nt0err, nfitmom );
    }

    if  ( (fitstatus == 1) ){   //fitstatus should be 1 for all the acceptable tracks

      bool inBins = (fitmom >= lowestBin && fitmom < topOfLastBin);
      bool minQcuts = false;
      
      if  ( (nactive   >= minimum_nactive) && (t0err     <= maximum_t0err)
	    && (fitmomerr <= maximum_fitmomerr)  && (fitcon    >= minimum_fitcon) 
	    && (fitmom    < 120.0) )minQcuts=true; 
      
      
      if(apply_timeCuts){//apply_time cut
	if(rg.enabled){ //if rg enabled
	  const int nhist = rg.nhist;
	  if (dio) {
	    if (use_diowt) { rg.hDIOspectrum->Fill(fitmom, diowt);   if(inBins)rg.hDIOspectrumX->Fill(fitmom, diowt);
	      for(int jj=0; jj<nhist; ++jj)rg.h_fitmomDIO[jj]->Fill(fitmom, diowt);
	      if(minQcuts)for(int jj=0; jj<nhist; ++jj)rg.hQ_fitmomDIO[jj]->Fill(fitmom, diowt); //fill hist with minQ cuts
	    } 
	    else {           rg.hDIOspectrum->Fill(fitmom, 1.0);	  rg.hDIOspectrumX->Fill(fitmom, 1.0);  
	      for(int jj=0; jj<nhist; ++jj)rg.h_fitmomDIO[jj]->Fill(fitmom, 1.0);
	      if(minQcuts)for(int jj=0; jj<nhist; ++jj)rg.hQ_fitmomDIO[jj]->Fill(fitmom, 1.0); //fill hist with minQ cuts
	    }
	  } //use dio
	  else {             rg.hCEspectrum->Fill(fitmom, 1.0);       rg.hCEspectrumX->Fill(fitmom, 1.0); 
	    for(int jj=0; jj<nhist; ++jj)rg.h_fitmomCE[jj]->Fill(fitmom, 1.0);
	    if(minQcuts)for(int jj=0; jj<nhist; ++jj)rg.hQ_fitmomCE[jj]->Fill(fitmom, 1.0); //fill hist with minQ cuts
	  }
	}
      }

      if(!minQcuts)continue;
      for (size_t t = 0; t < nt; ++t) {
	if ( (!apply_timeCuts) || (t0 >= tCuts[t]) ) {   // found a good track!
	  ++tracksExtracted;
	  double p = fitmom;
          if (inBins) { // p in range to be binned
            unsigned int binNumber = (p-lowestBin)/binSize;
            if (binNumber >= nBins) binNumber = nBins-1;
            if ( use_diowt ) {
              counts[t][binNumber] += diowt;  // weighted DIO tracks
              dioTotalWeight[t] += diowt;
              if ( fitmom >= canonicalRangeLo && fitmom <= canonicalRangeHi ) {
                dioTotalWeightInCanonicalRange[t] += diowt;
              }
            } else {
              counts[t][binNumber] += 1.0; // everything has equal weight
            }
            // #define RPCDIAGNOSTIC
#ifdef  RPCDIAGNOSTIC
            if ( fitmom >= canonicalRangeLo && fitmom <= canonicalRangeHi ) {
              if (!apply_timeCuts) {
                std::cout << p << "  " << t0 << "\n";
              }
            }                 
#endif           
            ++nBinnedTracks[t];
            ++tracksBinned;
            if (dio) highPs[t].add(p);
          } // end of *if* on p in range to be binned
        } // end of *if* on found a good track
      } // end of *for* on t going up to nt 
    } // end of *if* checking the track quality cuts    
  } // end of *for* loop on i which goes thru track entries
  //std::cerr << "### 8\n";
  if (OUTPUT_fileEntriesStatistics) {
    os << "There are " << ntracks << " entries in the tree\n";
    for ( size_t t = 0; t < nt; ++t ) {
      os << "There are " << nBinnedTracks[t] 
         << " binned entries passing the cuts with t0 >= " 
         << tCuts[t] << "\n";
      if (use_diowt) { 
        os << "(total of DIO weights passing cuts with t0 >= "
           <<  tCuts[t] << " is " << dioTotalWeight[t] << ")\n";
        os << "total weight in " << canonicalRangeLo << " - " 
           << canonicalRangeHi << " is " << dioTotalWeightInCanonicalRange[t]
           << "\n";
      }
    }
    if (OUTPUT_detailedCutTrace){
      outputCutSequence ( treeEntries, nfitstatus, nnactive, nfitcon, 
			  nfitmomerr, nt0err, nfitmom );
    }
  }
              
} // extractFitmom

void FMtool::extractFitmom 
  ( std::vector<std::string> const & listOfFileNames
  , std::vector< std::vector<double> > & counts
  , bool dio
  , bool use_diowt
  , bool apply_timeCuts 
  )
{
  // Obtain the momentum data as a vector of momenta
  if (OUTPUT_fileNames) {
    os << "extracting fitmom from a chain of files \n";
  }
  size_t nt = tCuts.size();
  if (dio) {
    for (size_t i = 0; i < nt; ++i) {
      NthHighest nl(minimumMeaningfulDIOtail, topOfLastBin);
      highPs.push_back(nl);
    }
  }
  for (unsigned int i=0; i<listOfFileNames.size(); ++i) {
    if (OUTPUT_fileNames) {
      os << "["  << i << "] : "
                << listOfFileNames[i].c_str() << "\n";
    }
    TTreeAccessor ta = accessTTree ( listOfFileNames[i] );
    //std::cerr << "### TTreeAccessor formed\n";
    extractFitmom 
      (  ta.tracks,  counts, highPs, dio, use_diowt, apply_timeCuts );
    //std::cerr << "### extractFitmom(ta.tracks...) returned \n";
    delete ta.tracks;   ta.tracks = 0;
    delete ta.tfp;      ta.tfp = 0;
  }
  if (dio) {
    lowMomentumCutCeiling.clear();
    lowMomentumCutCeiling.resize(nt);
    for ( size_t t = 0; t < nt; ++t ) {
      lowMomentumCutCeiling[t] = highPs[t].value();
      if (OUTPUT_fileEntriesStatistics) {
        os << minimumMeaningfulDIOtail 
           << "-th highest momentum of binned good tracks with t0 >= "
           << tCuts[t] << " is " 
           << lowMomentumCutCeiling[t] << "\n";
      }
    }
  }
} 

FMtool::TTreeAccessor 
FMtool::accessTTree(std::string const & fileName) const
{
  TTreeAccessor ta;
  { ifstream test(fileName.c_str());
    if (!test) {
      std::cerr << "Cannot open file " << fileName.c_str() << "\n";
      std::exit(1);
    }
  } 
  ta.tfp  = new TFile (fileName.c_str());

  if (ta.tfp->IsZombie()) {
    std::cerr << "Failed to properly open the output file " 
              << fileName << "\n";
    std::exit(1);
  }
  //  std::cout << "tft is not a zombie \n";
  ta.tracks = (TTree*)ta.tfp->Get ("ReadKalFits/trkdiag");
  if (ta.tracks != 0) {
    if (OUTPUT_tracksTTreeLocation) 
        os << "(tracks TTree found in ReadKalFits/trkdiag)\n";
  } else {
    ta.tracks = (TTree*)ta.tfp->Get ("trkdiag");
    if (ta.tracks != 0) {
      if (OUTPUT_tracksTTreeLocation) 
        os << "(tracks TTree found in trkdiag)\n";
    }
  }
  if (ta.tracks == 0) {
    std::cerr << "In file " << fileName << ":\n"
              << "  Failed to open tracks at "  
              << "ReadKalFits/trkdiag or trkdiag" <<  "\n";
    os        << "In file " << fileName << ":\n"
              << "  Failed to open tracks at "  
              << "ReadKalFits/trkdiag or trkdiag" <<  "\n";
    std::exit(1);
  }
  return ta;
} // accessTTree

void FMtool::accumulateCutStatistics( int fitstatus, int nactive,  double fitcon,
				      double fitmomerr, double t0err, double fitmom, 
				      int & nfitstatus,  int & nnactive, 
				      int & nfitcon, int & nfitmomerr, 
				      int & nt0err, int & nfitmom ) const 
{
  if (fitstatus == 1) {
    ++ nfitstatus;
    if (nactive   >= minimum_nactive) {
      ++nnactive;
      if (fitcon    >= minimum_fitcon)  {
        ++nfitcon;
        if (fitmomerr <= maximum_fitmomerr) {
         ++nfitmomerr;
          if (t0err     <= maximum_t0err) {
            ++nt0err;
            if (fitmom    < 120.0) {
              ++nfitmom;
            }
          }
        }
      }
    }
  }
}

void FMtool::outputCutSequence ( int treeEntries, int nfitstatus, int nnactive,
		  int nfitcon, int nfitmomerr, int nt0err, int nfitmom ) const 
{
      os << "Cuts sequence: " << treeEntries << " tracks\n"
         << "passing fitstatus cut = " << nfitstatus << "  " 
         << (double)nfitstatus/treeEntries << "\n"
         << "passing nactive cut   = " << nnactive << "  " 
         << (double)nnactive/nfitstatus << "\n"
         << "passing fitcon cut    = " << nfitcon << "  " 
         << (double)nfitcon/nnactive << "\n"
         << "passing fitmomerr cut = " << nfitmomerr << "  " 
         << (double)nfitmomerr/nfitcon << "\n"
         << "passing t0err cut     = " << nt0err << "  " 
         << (double)nt0err/nfitmomerr << "\n"
         << "passing fitmom sanity = " << nfitmom << "  " 
         << (double)nfitmom/nt0err 
         << " ==> " <<  (double)nfitmom/treeEntries << "\n";
}

void FMtool::countMcmom 
  ( TTree * tracks, double & count103, double & count104, bool use_diowt )
{
  tracks->SetBranchStyle(0);
  tracks->SetBranchStatus("*",0);

  float mcmom; 
  if(tracks->GetBranch("mcmom")){
    tracks->SetBranchStatus("mcmom",1);  
    tracks->SetBranchAddress("mcmom" , &mcmom); 
  }  
  else {
    tracks->SetBranchStatus("mcinfo",1);
    //MCTrkInfo mcinfo;   
    tracks->SetBranchAddress("mcinfo", &mcinfo);
  }
  
  float diowt;  tracks->SetBranchStatus("diowt",1);  tracks->SetBranchAddress("diowt", &diowt);

  int treeEntries = tracks->GetEntries();
  CzarneckiDIOspectrumFunction cz(1.0);
  for (int i = 0; i < treeEntries; ++i) {
    int bytesRead = tracks->GetEntry(i);
    if(tracks->GetBranch("mcinfo"))mcmom=mcinfo._mom;
    if (bytesRead == 0) break;
    if  (  (mcmom >= 102.5) && (mcmom < 103.5) ) {
        count103 += 1.0;    
    } 
    if  (  (mcmom >= 103.5) && (mcmom < 104.5) ) {
        count104 += 1.0;    
    } 
    if (use_diowt) {
      double czweight = cz(mcmom);
      double wratio = diowt / czweight;
      if ( wratio < .99 || wratio > 1.01 ) {
        std::cerr << "Discrepancy in diowt: mcmom = " << mcmom
                  << " diowt = " << diowt 
                  << " expected " << czweight << "\n";  
      }
    }
  } 
} // countMcmom

void FMtool::countMcmom 
  ( std::vector<std::string> const & listOfFileNames
  , double & count103, double & count104, bool use_diowt )
{
  count103 = 0.0;
  count104 = 0.0;
  // Obtain the momentum data as a vector of momenta
  for (unsigned int i=0; i<listOfFileNames.size(); ++i) {
    TTreeAccessor ta = accessTTree ( listOfFileNames[i] );
    countMcmom 
      (  ta.tracks,  count103, count104, use_diowt );
    delete ta.tracks;   ta.tracks = 0;
    delete ta.tfp;      ta.tfp = 0;
  }
} 

void FMtool::normalizeDIObackground(double numberOfGeneratedDIOs, bool use_diowt) 
{
  // Normalization of the DIO background
  // Note that our convention is to normalize the background to a
  // number per proton on target, so the luminosity does not appear in 
  // the background normalization.  (But the luminosity is passed to the
  // FofM calculator because *it* needs a total number of background counts.)

  // All the DIO's in the files that generated DIO's according to the
  // czarnecki distribution were created from a distribution of DIO's 
  // with p >= 102.5 (or, in principle, some other lower limit).  If
  // the tracks were generated with a momentum distribution following
  // czarnecki, then the tracks should not be weghted, but an overall 
  // factor representing the fraction of DIO's above the minimum momnetum
  // must be used.  For a minimum p of 102.5, this factor is 4.42e-15.
  //
  // In other files, DIO's were generated flat in some momentum range 
  // (e.g., 95-105) in which case the tracks should be weighted according 
  // to diowt.  The typical diowt is in the area of 1.0e-15.

  double DIOgenerationWeighting = 1.0;
  
  if (use_diowt) {
    DIOgenerationWeighting = 
        DIOflatGenerationWindowHi - DIOflatGenerationWindowLo;
    if (OUTPUT_detailedTrace) os << "### DIOgenerationWeighting = "  << DIOgenerationWeighting << "\n";
  } else {
    DIOgenerationWeighting = fractionOfDIOsRepresented;
  }
 
  double DIOsPerStoppedMuon =  1.0 - capturedMuonsPerStoppedMuon;
  double DIObackgroundCountNormalization =
      ( stoppedMuonsPerPOT * DIOsPerStoppedMuon 
          * DIOgenerationWeighting  * liveGateFraction ) /
      ( numberOfGeneratedDIOs );

  if ( adHocDIOrescaler != 1.0 ) {
     os << "\n** For this trial, DIO normalization is scaled to "
        << adHocDIOrescaler << " times its actual value **\n\n";
     DIObackgroundCountNormalization *= adHocDIOrescaler;
  }

  os << "|| DIOgenerationWeighting G = " << DIOgenerationWeighting << "\n";
  if (use_diowt) {   
    os << "|| (" << DIOflatGenerationWindowHi << " - " << DIOflatGenerationWindowLo << ")\n";
  } else {
    os << "|| (fractionOfDIOsRepresented) \n";
  }
  os << "|| stoppedMuonsPerPOT     s = " << stoppedMuonsPerPOT << "\n";
  os << "|| DIOsPerStoppedMuon     d = " << DIOsPerStoppedMuon << "\n";
  os << "|| liveGateFraction       F = " << liveGateFraction << "\n";
  os << "||                  G*s*d*F = " << DIOgenerationWeighting*stoppedMuonsPerPOT*
                          DIOsPerStoppedMuon*liveGateFraction << "\n";
  os << "|| Divided by numberOfGeneratedDIOs (" << numberOfGeneratedDIOs
     << ") = " <<   DIOgenerationWeighting*stoppedMuonsPerPOT*
                          DIOsPerStoppedMuon*liveGateFraction 
                          / numberOfGeneratedDIOs << "\n"; 
  os << "|| DIObackgroundCountNormalization = " 
       << DIObackgroundCountNormalization << "\n";                   
  os << "|| Times Luminosity " << protonsOnTarget << " = " 
        << DIObackgroundCountNormalization * protonsOnTarget << "\n";             
  if (use_diowt) {   
    os << 
    "|| (This is number of DIO counts in experiment per unit DIO weight in input files) \n";
  } else {
    os << "|| (This is number of DIO counts per DIO in the input files) \n";
  }

  DIObackgroundNormalization = DIObackgroundCountNormalization/binSize;
  for (size_t t =0; t <  tCuts.size(); ++t) {
    os << "Normalizing DIO background for t = " << tCuts[t] << "\n";
    normalizeAbackground( "DIO", 
      DIObackgroundCountNormalization, DIObackground[t] );
  }

} // normalizeDIObackground


void FMtool::normalizeRPCbackground(double numberOfGeneratedRPCs) 
{
  // Normalization of the RPCbackground
  // Note that our convention is to normalize the background to a
  // number per proton on target, so the luminosity does not appear in 
  // the background normalization.

  // Explanation of factors considered:
  //
  // We start by forming RPCnormalization which is the number of RPC
  // electrons over the lifetime of the experiment, represented by one
  // count in the files.  We illustrate the caluclation with numbers
  // taken from Bob's first RPC files.
  //
  // The normalization function has been passed  
  // (2.0e7).  So the files represent 2.0e7 "incidents," where an
  // incident is the forcing of a stopped pion to create one RPC electron.  
  // So one count represents 1 RPC per numberOfGeneratedRPCs incidents.
  //
  // Then we know that in reality, not every stopped pi leads to an RCP:
  // we have RPCperStoppedPion (.021) incidents per per stopped pion.
  // Thus there are RPCperStoppedPion/numberOfGeneratedRPCs RPC's per
  // stopped pion.
  //
  // Next, we have RPCtimeProbability stopped pions (**in the accepted time
  // window**) per POT.  For example, with t0 = 792, RPCtimeProbability was
  // 2.5e-16.  This probability is a combination of the chance of generating a 
  // pion that stops in the target, times the chance of that pion surviving
  // to at least t0.
  // Thus (although we will multiply by RPCtimeProbability later on) we have
  // RPCperStoppedPion/numberOfGeneratedRPCs * RPCtimeProbability
  // RPC electrons per POT.
  // 
  // Since we want a number per POT, we need not factor in 
  // the number of protons on target.  RPCbackgroundCountNormalization is
  //
  // 1 count = 1 RCP / numberOfGeneratedRPCs 
  //         * RPCperStoppedPion * RPCtimeProbability 
  //
  // If we want, we can normalize RPCnormalization to POT times this,
  // saying how many RPC's each count represents over the experiment lifetime.
  // With the original RPC files this was 1/(2e7) * .021 * 2.5e-16 * 3.6e20
  // which was 9.45e-5.

  // The number 9.45e-5 was based on just the processes Bob B. had
  // simulated to get these files.  But in fact there are other processes,
  // in particular, internal conversion, which was left out is is about
  // as intense as RPC photon conversions in the stopping target.
  // So at the time we multiplied by a further nternalRPCconversionFactor of 2.
  // However, we will assume now that all the relevant processes are included
  // in the simulation.  (If not, one can supply an ad hoc factor for RPC's.)

  //ParameterSet RPCpset = 
  //  pset.get<ParameterSet>("radiativePionCaptureElectrons");  

  double RPCnormalizationBase =  RPCperStoppedPion / numberOfGeneratedRPCs;  
  // RPCnormalizationBase still needs a factor of RPCtimeProbability

  if ( adHocRPCrescaler != 1.0 ) {
     os << "\n**** For this trial, RPC normalization is scaled to "
              << adHocRPCrescaler << " times its actual value ****\n\n";
     RPCnormalizationBase *= adHocRPCrescaler;
  }

  for (size_t t = 0; t <  tCuts.size(); ++t) {
    double RPClivePionFraction =  RPCtimeProbability(tCuts[t]);
    os << "|| stoppedPionsPerPOT                    s = " << stoppedPionsPerPOT << "\n";
    os << "|| fraction of pions lasting to t = " << tCuts[t] << ": F = " 
       << RPCtimeProbability(tCuts[t])/stoppedPionsPerPOT << "\n";
    os << "|| surviving stopped pions from proton pulse           s*f = "  
       << RPCtimeProbability(tCuts[t]) << "\n";
    double flatPionBackgroundInLiveGate = extinction * stoppedPionsPerPOT 
                * (cadence - tCuts[t]) / cadence;
        // This is the number of stopped pions arriving in the live gate
        // (between the cut time and the end of the time window) per
        // desirable POT, due to the (tiny) number of POT's outside the
        // proton pulse (the non-extinct protons)   
    RPClivePionFraction += flatPionBackgroundInLiveGate;
    os << "|| stoppedPionsPerPOT                    s = " << stoppedPionsPerPOT << "\n";
    os << "|| extinction factor                     e = " << extinction << "\n";
    os << "|| live RPC gate fraction past t = " << tCuts[t] << ":  g = "  
       << (cadence - tCuts[t]) / cadence << "\n";
    os << "|| surviving stopped pions from flat background      s*e*g = " 
       <<  flatPionBackgroundInLiveGate << "\n";
    os << "|| sum of stopped pions past t = "  << tCuts[t] << ":      T = s*f + s*e*g = "
       <<  RPCtimeProbability(tCuts[t]) + flatPionBackgroundInLiveGate << "\n";
    os << "|| RPCperStoppedPion     r = " << RPCperStoppedPion << "\n";
    double RPCperPOTpastT = RPCperStoppedPion * RPClivePionFraction;
    os << "|| RPCs per POT past t = "  << tCuts[t] 
       << ":                r*T =  " <<  RPCperPOTpastT << "\n";
    if ( adHocRPCrescaler != 1.0 ) {
      RPCperPOTpastT *= adHocRPCrescaler;       
      os << "|| Times an ad hoc RPC rescaling factor of " << adHocRPCrescaler 
          << " --> r*T = " << RPCperPOTpastT << "\n";
    }
    os << "|| Divided by numberOfGeneratedRPCs (" << numberOfGeneratedRPCs
       << ") = " <<  RPCperPOTpastT / numberOfGeneratedRPCs << "\n"; 
    double RPCbackgroundCountNormalization = 
                RPCnormalizationBase * RPClivePionFraction;
    os << "|| RPCbackgroundCountNormalization = " << RPCbackgroundCountNormalization << "\n";     os << "|| Times Luminosity " << protonsOnTarget << " = " 
        << RPCbackgroundCountNormalization * protonsOnTarget << "\n";             
    os << "|| (This is number of RPC counts in the experiment per RPC in the input files) \n";

    if (OUTPUT_RPClivePionFraction) {
      os << "RPClivePionFraction = " << RPClivePionFraction << " + " 
         << flatPionBackgroundInLiveGate << " = " 
         << RPClivePionFraction + flatPionBackgroundInLiveGate << "\n";
    }
    
    if (OUTPUT_backgroundStrength) {
      os << "\nWith minimum t0 = " << tCuts[t]
         << " and extinction = " << extinction << "\n"
         << "RPCliveGateFraction is "  
         << RPClivePionFraction - extinction * stoppedPionsPerPOT
         << " + " << extinction * stoppedPionsPerPOT 
         << " = " << RPClivePionFraction << "\n"
         << "RPCbackgroundCountNormalization is " 
         << RPCbackgroundCountNormalization << "\n";
    }

    normalizeAbackground( "RPC", 
      RPCbackgroundCountNormalization, RPCbackground[t] );
  }

} // normalizeRPCbackground

void FMtool::normalizeAbackground(  
        std::string const & bname, 
        double backgroundCountNormalization, 
        std::vector<double> & background)
{        
  double backgroundNormalization =  backgroundCountNormalization / binSize;  
             
  if (OUTPUT_backgroundStrength) {
    os <<  bname << "background: (Count normalization is " 
       << backgroundCountNormalization << " ==> " 
       << backgroundCountNormalization*protonsOnTarget << ")\n";
  }

  double canonicalRangeCount = 0;
  for (unsigned int i=0; i < background.size(); ++i) {
    double p = lowestBin + binSize*i;
    if ( p >= canonicalRangeLo && p < canonicalRangeHi )  {
      canonicalRangeCount += background[i];
    }
  }
  if (OUTPUT_backgroundStrength) {
    os << "There are " <<  canonicalRangeCount
       << " " << bname << "'s in the files in the range of "
       << canonicalRangeLo << " - "  << canonicalRangeHi << " \n";
    os << "With normalization, this amounts to "
       << canonicalRangeCount * backgroundCountNormalization 
                              * protonsOnTarget 
       << " " << bname << " events in the range \n";
  }
  for (unsigned int i=0; i < background.size(); ++i) {
    background[i] *= backgroundNormalization;
    if (OUTPUT_backgroundStrength && OUTPUT_spectra) {
      os << lowestBin + binSize*i << ":   "
         << background[i] << " --> " 
         << background[i] / backgroundCountNormalization;
      os << " --> " << background[i] 
         << " ==> " << background[i] * protonsOnTarget << "\n";
    }
  }

} // normalizeAbackground

void FMtool::NthHighest::add(double p) 
{
  // Invariants:  vals is of size n; vals contains only supplied 
  // values of p and values top; for all m<k vals[m] <= vals[k];
  // m values have been supplied thus far.  Thus vals[0] is always
  // either the n-th highest or the highest if fewer than m were supplied.
  if (n==0) return;
  assert (p <= top);
  if (m==0) { 
    vals[0] = p;
  } else if (m < n) {
    double q = p;
    for (int i = 0; i<m; ++i) {
      if (q < vals[i]) {
        double t = vals[i];
        vals[i] = q;
        q = t;
      }
    } 
    vals[m] = q; 
  } else {
    if (p <= vals[0]) {
      ++m;
      return;
    }
    vals[0] = p;
    for (int i = 1; i<n; ++i) {
      if ( p > vals[i]) {
        vals[i-1] = vals[i];
        vals[i]   = p;
      }  else {
        break;
      }
    }
  }
  ++m;    
#ifdef DIAGNOSTIC_FOR_NTHHIGHEST
  if ( m < 40 && (p <= vals[0]  || m <= n) ) {
    std::cout << "~~~~~ p: " << p << "\n";
    for (int j = 0; j < n; ++j) {
      std::cout << "~~~~~ " << vals[j] << "\n";
    }
    std::cout << "\n\n";
  }
#endif
}

int FMtool::NthHighest::howManyHigher(double p) const 
{
  if (n== 0) return 0;
  // The following is inefficient by a factor of at least 2 but that is moot in this context 
  double k = 0;
  size_t vs = vals.size();
  for (size_t i = 0; i < vs; ++i) {
    if (vals[i] > p) ++k;
  }
  return k;
}

void FMtool::makeDIOreferenceTF1(int tracksExtracted, double DIOflatGenerationWindowLo)
{
  if (!rg.enabled) return;
  CzarneckiDIOspectrumFunction cz_unit(1.0);
  double czIntegral = cz_unit.integrated(DIOflatGenerationWindowLo);
  double normalization = tracksExtracted/czIntegral;
  CzarneckiDIOspectrumFunction cz(normalization);
  rg.fDIOreference = new TF1("fDIOreference",cz,DIOflatGenerationWindowLo,104.96, 0,
			     "Czarnecki-Tormo DIO spectrum");
  rg.rootfp->Append(rg.fDIOreference);
  rg.cint <<  "// Established TF1 function fDIOreference -- Czarnecki-Tormo DIO spectrum\n";
}


void FMtool::applyFofM() 
{
  //
  // Use the FofM class to calculate quantities of interest,
  // once for each t cut.
  //

  // Since our vectors represent values at the middle of intervals,
  // we should state our sample points vectors as such when using
  // them to construct spectra
  double midBinOffset = (topOfLastBin - lowestBin)/(2.0*nBins);
  double lowestPoint = lowestBin + midBinOffset;
  double endingPoint = topOfLastBin - midBinOffset;
  double midSplineBinOffset = (topOfLastBin - lowestBin)/(2.0*nSplineBins);
  double lowestSplinePoint = lowestBin + midSplineBinOffset;
  double endingSplinePoint = topOfLastBin - midSplineBinOffset;

  os << "\n------------------------------------ \n\n"
     << "Using the FofM tool \n\n";

  os << "lowestBin = " << lowestBin
     << "  topOfLastBin = " << topOfLastBin 
     << "\nlowestPoint = " << lowestPoint
     << "  endingPoint = " << endingPoint
     << "  nBins = " << nBins  << "\n"   
     << "nSplineBins = " << nSplineBins 
     << "  lowestSplinePoint = " << lowestSplinePoint
     << "  endingSplinePoint = " << endingSplinePoint << " \n\n";    

  // TODO - think out and implement the way we want users to exercise this
  //        control
  double lowCut;
  double highCut;
  ParameterSet pCutpset = pset.get<ParameterSet>("momentumCuts");
  bool optimize_pcut_low = pCutpset.get<bool>("optimize_pcut_low",true);
  lowCut = optimize_pcut_low ? 0 : pCutpset.get<double>("fixed_pcut_low");
  bool optimize_pcut_high = pCutpset.get<bool>("optimize_pcut_high",true);
  highCut = optimize_pcut_high ? 0 : pCutpset.get<double>("fixed_pcut_high");
        recordFclGroup("momentumCuts");
        recordFcl("optimize_pcut_low ",optimize_pcut_low);
	if (!optimize_pcut_low) recordFcl("fixed_pcut_low   ",lowCut);
        recordFcl("optimize_pcut_high",optimize_pcut_high);
	if (!optimize_pcut_high) recordFcl("fixed_pcut_high  ",highCut);
  std::string mfcString;
  MeritFunctionChoice mfc;
  // Obsolete, but here so that older fcl files will still work the same way: 
  bool useSmoothedPunziFMtoOptimize = pset.get<bool>
                        ("useSmoothedPunziFMtoOptimize", false);
  std::string mfcStringDefault = useSmoothedPunziFMtoOptimize ?
				  "SmoothedPunzi" : "FCsensitivity";
  // Newer way to specify merit function choice:
  mfcString = pset.get<std::string>("meritFunction","mfcStringDefault");
  if ( mfcString == "SmoothedPunzi" ) mfc = SmoothedPunziMeritFunction;
  if ( mfcString == "FCsensitivity" ) mfc = FCsensitivityMeritFunction;
        recordFclGroup("Optimization controls (not a separate parameter set)");
        recordFcl (  "meritFunction", mfcString );

  std::string table;
  std::vector<FofM::Summary> summaries(tCuts.size());
  double best_tCutMerit  = 0;
  size_t best_tCutNumber = 0;
  for (size_t tCutNumber=0; tCutNumber<tCuts.size(); ++tCutNumber) {
    //std::cerr << "##### summaries[" << tCutNumber << "]\n";
    summaries[tCutNumber] = 
        applyFofM( tCutNumber, lowestSplinePoint, endingSplinePoint, table, 
                        mfc, lowCut, highCut );
    if (summaries[tCutNumber].figureOfMerit > best_tCutMerit) {
      best_tCutMerit  = summaries[tCutNumber].figureOfMerit ;
      best_tCutNumber = tCutNumber;
    }
    if (OUTPUT_allTables) {
      os << table << "\n--------------------------------------- \n\n";
    }
  }

  os << "\nFigure of merit versus minimum of t0 window: \n\n";
  os << "tCut    FofM            SES         CL90         SPunzi    "
     << " pLow    Phigh \n";
  for (size_t tCutNumber=0; tCutNumber<tCuts.size(); ++tCutNumber) {
    os << tCuts[tCutNumber] 
       << "  " << summaries[tCutNumber].figureOfMerit;
    if (tCutNumber == best_tCutNumber) {
      os << " ** ";
      bottomLine("time cut",  tCuts[tCutNumber]); 
      bottomLine("lower momentum cut", summaries[tCutNumber].pCutLo);     
      bottomLine("upper momentum cut", summaries[tCutNumber].pCutHi);     
      bottomLine("Feldman/Cousins 90% sensitivity", summaries[tCutNumber].CL90sensitivity);   
      if(rg.enabled){
	int ivalv = 0;
	rg.valv[ivalv] = tCuts[tCutNumber];              ++ivalv;
	rg.valv[ivalv] = summaries[tCutNumber].pCutLo;   ++ivalv;
	rg.valv[ivalv] = summaries[tCutNumber].pCutHi;   ++ivalv;
	rg.valv[ivalv] = summaries[tCutNumber].CL90sensitivity;   ++ivalv;
	rg.valv[ivalv] = minimumMeaningfulDIOtail;   ++ivalv;
	rg.valv[ivalv] = canonicalRangeLo;   ++ivalv;
	rg.valv[ivalv] = canonicalRangeHi;   ++ivalv;
	rg.valv[ivalv] = lowestBin;   ++ivalv;
	rg.valv[ivalv] = topOfLastBin;   ++ivalv;
	rg.valv[ivalv] = stoppedMuonsDef;   ++ivalv;
	rg.valv[ivalv] = stoppedMuonsThisRun;   ++ivalv;
	rg.valv[ivalv] = sc_factor;   ++ivalv;
	rg.valv[ivalv] = stoppedMuonsPerPOTold;   ++ivalv;
	rg.valv[ivalv] = stoppedMuonsPerPOT;   ++ivalv;
	rg.tree->Fill();
      }
    } else {
      os << "    ";
    }
    os <<  summaries[tCutNumber].singleEventSensitivity 
       << "  " << summaries[tCutNumber].CL90sensitivity 
       << "  " << summaries[tCutNumber].smoothedPunziSensitivity;
    int pr = os.precision(5); 
    os << "  "  <<  summaries[tCutNumber].pCutLo
       << "  "  <<  summaries[tCutNumber].pCutHi << "\n\n";
    os.precision(pr);
  }

// Repeat the best one if there were more than 2 t cuts explored
// or if we have not been outputting tables:

#ifdef REMOVED_DIAGNOSTICS
  os << "tCuts.size() = " << tCuts.size() << " and OUTPUT_allTables = "
     << OUTPUT_allTables << "\n";
  summaries[best_tCutNumber] = 
        applyFofM( best_tCutNumber, lowestSplinePoint, endingSplinePoint, table,
                        mfc, lowCut, highCut );
        os << table << "\n ***** \n";
#endif
          
  if ( (tCuts.size() > 2) || (!OUTPUT_allTables) ) {
    //std::cerr << "***** summaries[" << best_tCutNumber << "]\n";
    summaries[best_tCutNumber] = 
        applyFofM( best_tCutNumber, lowestSplinePoint, endingSplinePoint, table, 
                        mfc, lowCut, highCut );
    //std::cerr << "***** summaries[" << best_tCutNumber << "] returned \n";
    os << table << "\n--------------------------------------- \n\n";
  }
  
  os << botMenu.str() << "\n";
  os << botLine.str() << "\n";

  os << "\n--------------------------------------- \n\n";
  
}  // applyFofM()

FofM::Summary
FMtool::applyFofM( size_t tCutNumber,
                   double lowestPoint,
                   double endingPoint,
                   std::string & table, 
                   MeritFunctionChoice mfc,
                   double lowCut, double highCut ) const
{
  //
  // Use the FofM class to calculate quantities of interest
  //


  os << "================================= \n"
     << "Time window minimum = " << tCuts[tCutNumber] << "\n"
     << "================================= \n";

  std::vector<double> sig = rebin(sigEfficiency[tCutNumber],nSplineBins,
                                lowestPoint, endingPoint);
  std::vector<double> dio = rebin(DIObackground[tCutNumber],nSplineBins,
                                lowestPoint, endingPoint);
  FofM::Spectrum CE_spectrum  (sig, lowestPoint, endingPoint);
  FofM::Spectrum DIO_spectrum (dio, lowestPoint, endingPoint);
  FofM figureOfMeritCalculator ( CE_spectrum, 
                                 DIO_spectrum, 
                                 protonsOnTarget,
                                 mfc, os );  

  if (OUTPUT_backgroundSplines) { 
    figureOfMeritCalculator.displayBackground(lowestBin, topOfLastBin, nBins);
    splines::Spline<1> sbkg = figureOfMeritCalculator.getBackground();
    splines::Grid<1> sgrid = figureOfMeritCalculator.getGrid();
    os << " L * DIO Background (normalized) vs L * Background (spline)\n";
    os << " The normalized numbers are counts in bin, per experiment\n";
    double roughInt = 0.0;
    double gridstep = sgrid[1] - sgrid[0];
//    double DIObackgroundCountNormalization = 
//      DIObackgroundNormalization * gridstep;
    for (unsigned int i = 0; i < sgrid.nPoints(); ++i) {
      double p =  sgrid[i];
      os << i << ": p = " << p  << "  " 
         << DIObackground[tCutNumber][i] / DIObackgroundNormalization 
         << "  "
         << protonsOnTarget * DIObackground[tCutNumber][i] << "  "
         << protonsOnTarget * sbkg(p) << "\n";
      if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) )  {
        roughInt += sbkg(p);
      }
    } 
    os << "Rough integral in canonical range based on sbkg = " 
       << protonsOnTarget*gridstep*roughInt << "\n";
  }

  // Add the RPC background - note that though time cuts are handled
  // analytically for RPC's, we are keeping one RPCbackground spectrum 
  // for each time cut.
  std::vector<double> rpc = rebin(RPCbackground[tCutNumber],nSplineBins,
                                lowestPoint, endingPoint);
  FofM::Spectrum RPC_spectrum (rpc, lowestPoint, endingPoint);
  figureOfMeritCalculator.addBackground (RPC_spectrum);

  if (OUTPUT_backgroundSplines) { 
    figureOfMeritCalculator.displayBackground(lowestBin, topOfLastBin, 100);
    splines::Spline<1> sbkg = figureOfMeritCalculator.getBackground();
    splines::Grid<1> sgrid = figureOfMeritCalculator.getGrid();
    os << " L * (DIO+RPC) Background (normalized) vs L * Background \n";
    double roughInt = 0.0;
    double gridstep = sgrid[1] - sgrid[0];
    for (unsigned int i = 0; i < sgrid.nPoints(); ++i) {
      double p =  sgrid[i];
      os << i << ": p = " << p  << "  " 
         << DIObackground[tCutNumber][i] / DIObackgroundNormalization 
         +  RPCbackground[tCutNumber][i] / RPCbackgroundNormalization << "  "
         << protonsOnTarget * 
           (DIObackground[tCutNumber][i] + RPCbackground[tCutNumber][i]) << "  "
         << protonsOnTarget * sbkg(p) << "\n";
      if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) )  {
        roughInt += sbkg(p);
      }
    } 
    os << "DIO + RPC Rough integral in canonical range based on sbkg = " 
              << protonsOnTarget*gridstep*roughInt << "\n";
  }

#ifdef EXERCISE_RESULTS_METHOD
  double sigEff;
  double bkgd;
  double figOfM;
  int    discoveryThresholdN;
  double discoveryThreshold;
  double smoothedWorstCaseSensitivity;
  double punziSensitivity;
  double CL90Sensitivity;
  double singleEventSensitivity;
  std::vector<double> upperLimitBR;
  std::vector<double> CLlowerBR;
  std::vector<double> CLupperBR;
    figureOfMeritCalculator.results 
      (maximumSignalCount, lowCut, highCut, sigEff, bkgd, figOfM,  
       discoveryThresholdN, discoveryThreshold, 
       smoothedWorstCaseSensitivity,punziSensitivity,
       CL90Sensitivity, singleEventSensitivity,
       upperLimitBR, CLlowerBR, CLupperBR); 
#endif

  FofM::Summary summary;
  std::string t; 
  if ((lowCut == 0) && (highCut == 0)) {    
    t = figureOfMeritCalculator.tables
                        (maximumSignalCount, lowCut, highCut, os, summary);
    if (lowCut > lowMomentumCutCeiling[tCutNumber]) {
      DIOstatisticsWarning(lowCut,tCutNumber);
      lowCut = lowMomentumCutCeiling[tCutNumber];
      t = figureOfMeritCalculator.tables_fixed_lowCut 
                        (maximumSignalCount, lowCut, highCut, os, summary);
    }
  } else if (lowCut == 0) { 
    t = figureOfMeritCalculator.tables_fixed_highCut 
                        (maximumSignalCount, lowCut, highCut, os, summary);
    if (lowCut > lowMomentumCutCeiling[tCutNumber]) {
      DIOstatisticsWarning(lowCut,tCutNumber);
      lowCut = lowMomentumCutCeiling[tCutNumber];
      t = figureOfMeritCalculator.tables_fixed_cuts 
                        (maximumSignalCount, lowCut, highCut, os, summary);
    }
  } else if (highCut == 0) { 
    t = figureOfMeritCalculator.tables_fixed_lowCut 
                        (maximumSignalCount, lowCut, highCut, os, summary);
    if (lowCut > lowMomentumCutCeiling[tCutNumber]) {
      DIOstatisticsWarningFixedLowCut(lowCut,tCutNumber);
    }
  } else {
    t = figureOfMeritCalculator.tables_fixed_cuts 
                        (maximumSignalCount, lowCut, highCut, os, summary);
    if (lowCut > lowMomentumCutCeiling[tCutNumber]) {
      DIOstatisticsWarningFixedLowCut(lowCut,tCutNumber);
    }
  }

  splines::Spline<1> ceSpline  = CE_spectrum.representation();
  splines::Spline<1> dioSpline = DIO_spectrum.representation();
  splines::Spline<1> rpcSpline = RPC_spectrum.representation();
  
  double ceN  = ceSpline.integrate(lowCut, highCut) * protonsOnTarget * 1.0e-16;
  double dioN = dioSpline.integrate(lowCut, highCut) * protonsOnTarget;
  double rpcN = rpcSpline.integrate(lowCut, highCut) * protonsOnTarget;

  if (dioN < 0) {
    double splineBinSize = (topOfLastBin-lowestBin)/nSplineBins;
    double keyBinLeftSide = std::floor((lowCut - lowestBin)/splineBinSize); 
    // double midSplineBinOffset = (topOfLastBin - lowestBin)/(2.0*nSplineBins);
    // Perhaps I should be  adding midSplineBinOffset.  But this works.
    double newFixedLowCut = lowestBin + keyBinLeftSide*splineBinSize;
    DIOnegativeExplanation ( lowCut, tCutNumber, newFixedLowCut );
    t = figureOfMeritCalculator.tables_fixed_cuts 
                        (maximumSignalCount, newFixedLowCut, highCut, os, summary);
    lowCut = newFixedLowCut;
    ceN  = ceSpline.integrate(lowCut, highCut) * protonsOnTarget * 1.0e-16;
    dioN = dioSpline.integrate(lowCut, highCut) * protonsOnTarget;
    rpcN = rpcSpline.integrate(lowCut, highCut) * protonsOnTarget;
  }
  os << "\n In momentum range " << lowCut << " -- " << highCut << " There are: \n"
     << dioN << " DIO \n"
     << rpcN << " RPC \n"
     << ceN  << " CE at a branching ratio of 1.0e-16 \n";
  table = t;
   
  return summary;
  
} // applyFofM(cutNumber...)

std::vector<double> FMtool::rebin(std::vector<double> const & x, unsigned int n, 
				  double a, double b)
{
  bool diag = false;
  unsigned int m = x.size();
  if (m <= n) return x;
  std::vector<double> y(n);
  double s = a; // start of the active bin in loop over original binnning
  double t = a; // start of the active bin in loop over new binning
                // invariant:  s >= t at all times
  double ds = (b-a)/m;
  double dt = (b-a)/n;
  double sum = 0; 
  unsigned int j = 0;
  double rescale = (static_cast<double>(n))/m;
  if (diag) std::cout << "\n  rebin -- "  << m << " --> " << n << ": " << rescale << " \n\n";
  for (unsigned int i = 0; i != m; ++i) {
    if ( s + ds >= t + dt ) {
      sum += x[i] * (t + dt -s) / ds;
      y[j++] = sum*rescale;
      if (diag) std::cout << "a s = " << s << " t = " << t << " x[" << i << "] = " << x[i]
                          << " y[" << j-1 << "] " << y[j-1] << "\n";
      sum = 0;      
      if (j >= n) break;
      t += dt;      
      sum = x[i] * (s + ds -t) / ds;
      if (diag) std::cout << "b s = " << s << " t = " << t << " x[" << i << "] = " << x[i]
                          << " y[" << j << "] " << y[j] << "\n";
    } else {
      sum += x[i];
      if (diag) std::cout << "c s = " << s << " t = " << t << " x[" << i << "] = " << x[i]
                          << " sum = " << sum << "\n";
    }
    s += ds;
  }
  y[n-1] += sum*rescale;
  if (diag) std::cout << "d sum = " << sum << " y[" << n-1 << "] " << y[n-1] << "\n";
  return y;
}


void FMtool::DIOstatisticsWarning
                (double computedLowCut, size_t tCutNumber) const
{
  os << "Warning -- Low Statistics in DIO tail using time cut at " 
     << tCuts[tCutNumber] << "\n";
  os << "Explanation: \n"
     << "The automated momentum cut optimizer selected a lower p limit of "
     << computedLowCut << ".\n"
     << "However, there are " << highPs[tCutNumber].howManyHigher(computedLowCut) 
     << " < " << minimumMeaningfulDIOtail 
     << " DIO tracks that pass that momentum cut.\n"
     << "Sensitivity calculations based on so few tracks are deemed "
     << "untrustworthy.\n"
     << "The calculations will be redone fixing a lower momentum cut at "
     <<  lowMomentumCutCeiling[tCutNumber] << ".\n"
     << "The resulting FC90 and Smooth Punzi sensitivity numbers should be \n"
     << "regarded as conservative (high) values due to the low statistics.\n"
     << "The single event sensitivity should be regarded as optimistic \n" 
     << "due to the expanded momentum acceptance interval.\n\n";
}

void FMtool::DIOstatisticsWarningFixedLowCut
                (double fixedLowCut, size_t tCutNumber) const
{
  os << "Warning -- Low Statistics in DIO tail using time cut at " 
     << tCuts[tCutNumber] << "\n";
  os << "Explanation: \n"
     << "A lower p limit of " << fixedLowCut << " was specified.\n"
     << "There are " << highPs[tCutNumber].howManyHigher(fixedLowCut) 
     << " < " << minimumMeaningfulDIOtail 
     << " DIO tracks that pass that momentum cut.\n"
     << "Sensitivity calculations based on so few tracks are deemed "
     << "untrustworthy.\n\n";
}

void FMtool::DIOnegativeExplanation
(double computedLowCut, size_t tCutNumber, double newFixedLowCut) const
{
  os << "Note - the original automated optimization has suggested a lower p cut of "
     << computedLowCut << ".\n"
     << "  This is higher than the bin center for the highest momentum DIO in the data set.\n"
     << "  This situation be an artifact of eliminating too many of the actual causes of \n"
     << "  momentum uncertainty (for example, no proton absorber *and* no noise mix), \n"
     << "  and can result in a fictitious negative integral of the spline representing the \n"
     << "  the DIO part of the background.  \n"
     << "  We will deal with this by moving the low p cut down to a fixed cut at "
     << newFixedLowCut << "\n"
     << "  which is the center of the last bin containing non-zero DIO counts \n"
     << "  (the high cut stays as orignally determined). \n"; 
}


void FMtool::bottomLine (std::string const & description, double value) 
{
  botMenu << description << "\n";
  botLine << value << "  ";
} 

std::vector<double> FMtool::RPCtimeProbabilityVector() const {
  // returns the integrated (to t=1000) probability for 
  // pions for stopping times from 400 to 1000
  // For example, an entry of 9.3e-10 for t = 500 means 
  // that for a given proton on production target, the chance
  // of a pion stopping in the stopping target at a time AFTER
  // 500 nsec is 9.3 parts in 10^10.
  std::vector<double> v;
  v.push_back( 6.35462e-8 ); // 400 
  v.push_back( 4.49688e-8 );  
  v.push_back( 3.31549e-8 );  
  v.push_back( 2.13442e-8 );  
  v.push_back( 1.43251e-8 );  
  v.push_back( 9.45227e-9 );  
  v.push_back( 6.13425e-9 );  
  v.push_back( 3.91704e-9 );  
  v.push_back( 2.46267e-9 );  
  v.push_back( 1.52597e-9 );  
  v.push_back( 9.33404e-10 ); // 500 
  v.push_back( 5.64770e-10 );  
  v.push_back( 3.38856e-10 );  
  v.push_back( 2.02111e-10 );  
  v.push_back( 1.12106e-10 );  
  v.push_back( 7.12319e-11 );  
  v.push_back( 4.22078e-11 );  
  v.push_back( 2.50041e-11 );  
  v.push_back( 1.48127e-11 );  
  v.push_back( 8.77516e-12 );  
  v.push_back( 5.19851e-12 ); // 600 
  v.push_back( 3.07967e-12 );  
  v.push_back( 1.82446e-12 );  
  v.push_back( 1.08086e-12 );  
  v.push_back( 6.40343e-13 );  
  v.push_back( 3.79379e-13 );  
  v.push_back( 2.24782e-13 );  
  v.push_back( 1.33196e-13 );  
  v.push_back( 7.89386e-14 );  
  v.push_back( 4.67944e-14 );  
  v.push_back( 2.77499e-14 ); // 700 
  v.push_back( 1.64655e-14 );  
  v.push_back( 9.77835e-15 );  
  v.push_back( 5.81445e-15 );  
  v.push_back( 3.46343e-15 );  
  v.push_back( 2.06810e-15 );  
  v.push_back( 1.23909e-15 );  
  v.push_back( 7.45480e-16 );  
  v.push_back( 4.50540e-16 );  
  v.push_back( 2.73454e-16 );  
  v.push_back( 1.66537e-16 ); // 800 
  v.push_back( 1.01646e-16 );  
  v.push_back( 6.21046e-17 );  
  v.push_back( 3.79546e-17 );  
  v.push_back( 2.31921e-17 );  
  v.push_back( 1.41678e-17 );  
  v.push_back( 8.65286e-18 );  
  v.push_back( 5.28415e-18 );  
  v.push_back( 3.22705e-18 );  
  v.push_back( 1.97068e-18 );  
  v.push_back( 1.20319e-18 ); // 900 
  v.push_back( 7.34185e-19 );  
  v.push_back( 4.47457e-19 );  
  v.push_back( 2.71993e-19 );  
  v.push_back( 1.64380e-19 );  
  v.push_back( 9.82723e-20 );  
  v.push_back( 5.75389e-20 );  
  v.push_back( 3.22688e-20 );  
  v.push_back( 1.64364e-20 );  
  v.push_back( 6.40187e-21 );  
  v.push_back( 0.0 );         // 1000 
  return v;
} // RPCtimeShapeVector()

splines::Spline<1> FMtool::RPCtimeProbabilitySpline()  const
{
  double first_time = 400.0;
  double last_time = 1000.0;
  double grid_spacing = 10.0; // nsec
  unsigned int npoints = (last_time - first_time)/grid_spacing + 1;
  std::vector<double> probs = RPCtimeProbabilityVector();
  if ( npoints != probs.size() ) {
    std::cerr << "\n???? Discrepancy between RPCtimeProbabilityVector.size() = "
              << probs.size() << " and expected size " << npoints << "\n\n";  
  }
  splines::Grid<1> g ( first_time, last_time + grid_spacing, npoints );
  splines::Spline<1> s ( probs, g );
  return s;
}

void FMtool::decideVerbosity()
{
  ParameterSet v = pset.get<ParameterSet>("verbosity");
  OUTPUT_signalEfficiency = v.get<bool>("signalEfficiency", false);
  OUTPUT_backgroundStrength = v.get<bool>("backgroundStrength", false);
  OUTPUT_tracksTTreeLocation = v.get<bool>("tracksTTreeLocation", false);
  OUTPUT_progress = v.get<bool>("progress", false);
  OUTPUT_fileNames = v.get<bool>("fileNames", false);
  OUTPUT_dataProperties = v.get<bool>("dataProperties", false);
  OUTPUT_RPClivePionFraction = v.get<bool>("RPClivePionFraction",false);
  OUTPUT_spectra = v.get<bool>("spectra", false);
  OUTPUT_backgroundSplines = v.get<bool>("backgroundSplines", false);
  OUTPUT_allTables = v.get<bool>("allTables", false);
  OUTPUT_fileEntriesStatistics = v.get<bool>("fileEntriesStatistics",false);
  OUTPUT_detailedTrace = v.get<bool>("detailedTrace",false);
  OUTPUT_detailedCutTrace = v.get<bool>("detailedCutTrace",false);
}
 
} // end namespace mu2e

static ParameterSet obtainPset(std::string const & parametersFile) 
{
  // The following is part of the mantra seen in 
  // artexternals/fhiclcpp/v2_16_03/example
  // It may or may not be strictly necessary in terms of 
  // getting to the fhicl file desired:
  putenv(const_cast<char*>("FHICL_FILE_PATH=./test:."));
  cet::filepath_lookup policy("FHICL_FILE_PATH");
  // TODO - maybe use filepath_lookup_nonabsolute
  // cet::filepath_lookup_nonabsolute policy("FHICL_FILE_PATH");
  
  // Get the table from the file
  fhicl::intermediate_table tbl;
  std::cout << "The name of the fhicl file is " << parametersFile << "\n";
  fhicl::parse_document(parametersFile, policy, tbl);  
  
  ParameterSet p;
  // convert to ParameterSet
  fhicl::make_ParameterSet(tbl, p);
  
  return p;
}

int main (int argc, char*argv[]) {
  if ( argc != 2 ) {
    std::cout << "This application requires exactly one argument: \n"
              << "The name of the controlling parameter set (.fcl) file. \n";
    exit (-1);
  }
  ParameterSet pset = obtainPset( argv[1] );
  //std::cerr << "### 1\n";
  std::string outputName = pset.get<std::string>("outputFile","cout"); 
  if (outputName == "cout") {
    std::cout << "Parameter set file:  " << argv[1] << "\n";
    mu2e::FMtool f (pset, std::cout) ;
    f.analyze();
  } else if (outputName == "cerr") {
    std::cerr << "Parameter set file:  " << argv[1] << "\n";
    mu2e::FMtool f (pset, std::cerr) ;
    f.analyze();    
  } else {
    std::ofstream os (outputName.c_str());
    os << "FMtool output file:  " << outputName << "\n";
    os << "Parameter set file:  " << argv[1] << "\n";
    mu2e::FMtool f (pset, os) ;
    f.analyze();
  }
  return 0;
}
