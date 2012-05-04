// FMtool.cpp

#include "FigureOfMerit/inc/FMtool.h"

// C++ includes.
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <cstdlib>

// art and art externals includes

#include "fhiclcpp/ParameterSet.h"

// mu2e includes
#include "FigureOfMerit/inc/FofM.h"
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

namespace mu2e {

FMtool::FMtool() 
  : OUTPUT_signalEfficiency     (true)
  , OUTPUT_spectra              (false)
  , OUTPUT_backgroundStrength   (true)
  , OUTPUT_backgroundSplines    (false)
  , RPCtimeProbability(RPCtimeProbabilitySpline())
{
}

void FMtool::analyze()
{
  // TODO - replace the setup steps by a driver that uses a file or whatever
  //        to get data inputs.
  
  setFilterLevels();
  setCanonicals();
  double nCE = obtainCEdata();
  normalizeSignalEfficiency(nCE);  
  bool (use_diowt);
  double genDIO = obtainDIOdata(use_diowt);
  normalizeDIObackground(genDIO,use_diowt);  
  double genRPC = obtainRPCdata();
  normalizeRPCbackground(genRPC);  
  applyFofM();
}

void FMtool::setFilterLevels() 
{
  // TODO - get this data from a file or from the user; allow choices.

  TrackerCutsSet trackerCutsChoice = TrackerCuts_B;
  std::string choice;
  bool trackerCutsChoiceMade = false;
  while (!trackerCutsChoiceMade) {
    std::cout << 
      "Enter uppercase letter A, B, C, or D to specify tracker cuts choice: ";
    std::cin >> choice;
    trackerCutsChoiceMade = true;
    if        (choice[0] == 'A') {
      trackerCutsChoice = TrackerCuts_A;
    } else if (choice[0] == 'B') { 
      trackerCutsChoice = TrackerCuts_B;
    } else if (choice[0] == 'C') { 
      trackerCutsChoice = TrackerCuts_C;
    } else if (choice[0] == 'D') { 
      trackerCutsChoice = TrackerCuts_D;
    } else {
      std::cout << 
          "Tracker cuts choice must be A, B, C, or D\n";
      trackerCutsChoiceMade = false;
    }
  }
  ourTrackerCuts (trackerCutsChoice,       
                  minimum_nactive,
                  maximum_t0err,
                  maximum_fitmomerr,
                  minimum_fitcon);
  if ((trackerCutsChoice != TrackerCuts_D) && choice.size() > 1) {
    nextTrackerCuts(trackerCutsChoice);   
    double fraction = 0.0;
    if ( choice[1] >= '0' && choice[1] <= '9') { 
      fraction = (choice[1] - '0')/10.0;
    }
    int n_nactive = 0;
    float n_t0err = 0;
    float n_fitmomerr = 0;
    float n_fitcon = 0;
    ourTrackerCuts (trackerCutsChoice,       
                    n_nactive,
                    n_t0err,
                    n_fitmomerr,
                    n_fitcon);
    minimum_nactive = static_cast<int> 
        (fraction*n_nactive + (1.0-fraction)*minimum_nactive);
    maximum_t0err     = fraction*n_t0err     + (1.0-fraction)*maximum_t0err;
    maximum_fitmomerr = fraction*n_fitmomerr + (1.0-fraction)*maximum_fitmomerr;
    minimum_fitcon = std::exp
      ( fraction*std::log(n_fitcon) + (1.0-fraction)*std::log(minimum_fitcon) );
    std::cout << "minimum_nactive   = " << minimum_nactive << "\n"
              << "maximum_t0err     = " << maximum_t0err << "\n"
              << "maximum_fitmomerr = " << maximum_fitmomerr << "\n"
              << "minimum_fitcon    = " << minimum_fitcon << "\n";
  }

  std::cout << "Enter the t0 cut (or 0 to use pre-set live window fractions): ";
  std::cin >> minimum_t0;

}  // setFilterLevels

void FMtool::ourTrackerCuts ( TrackerCutsSet cuts, 
                              int & nactive, 
                              float & t0err, 
                              float & fitmomerr,
                              float & fitcon )
{
  switch (cuts) {
    case TrackerCuts_A:
      nactive   = 20;           // minimum
      t0err     = 10.0;         // maximum
      fitmomerr = 1.0;          // maximum
      fitcon    = 1.0e-6;       // minimum
    break;
    case TrackerCuts_B:
      nactive   = 20;           // minimum
      t0err     = 1.5;          // maximum
      fitmomerr = 0.2;          // maximum
      fitcon    = 1.0e-4;       // minimum
    break;
    case TrackerCuts_C:
      nactive   = 25;          // minimum
      t0err     = 1.0;         // maximum
      fitmomerr = 0.18;        // maximum
      fitcon    = 1.0e-3;      // minimum
    break;
    case TrackerCuts_D:
      nactive   = 30;          // minimum
      t0err     = 0.9;         // maximum
      fitmomerr = 0.15;        // maximum
      fitcon    = 1.0e-2;      // minimum
    break;
    default:
      std::cout << "???? trackerCutsChoice switch was not one of the cases!\n";
      std::cout << "???? using cuts set B by default\n";
      nactive   = 20;
      t0err     = 1.5; 
      fitmomerr = 0.2;
      fitcon    = 1.0e-4;
  }
} // ourTrackerCuts

void FMtool::nextTrackerCuts( TrackerCutsSet & cuts ) 
{
  TrackerCutsSet result;
  switch (cuts) {
    case TrackerCuts_A:
      result = TrackerCuts_B;
    break;
    case TrackerCuts_B:
      result = TrackerCuts_C;
    break;
    case TrackerCuts_C:
      result = TrackerCuts_D;
    break;
    default:
      result = cuts;
  }
  cuts = result;
}  // nextTrackerCuts

void FMtool::setCanonicals() 
{
  // TODO - allow this to be set from file or whatever.

  // TODO - some way to control which merit figure is used to optimize

  // Canonical range
  canonicalRangeLo = 103.51;
  canonicalRangeHi = 104.70;

  // binning info
  lowestBin = 102.5;  
  topOfLastBin = 106.0;  // top of last bin
  nBins = 350;
  binSize = (topOfLastBin - lowestBin)/nBins;
    
  // Experiment quantities
    // TODO -- obtain these normalization-related constants from some sort of 
    // database or conditions service. 
  protonsOnTarget = 3.6e20;  
  stoppedMuonsPerPOT = 0.0021;  // I have seen .00215; new  number is .0016
  capturedMuonsPerStoppedMuon = 0.609; // DocDB 48 - I have seen .59
  RPCperStoppedPion = 0.021;     // DocDB 1087
  stoppedPionsPerPOT = 1.53e-6;
  
  liveGateFraction = .47;  // or .50 -- captures in time window in 48
  if ( minimum_t0 != 0 )  liveGateFraction = 1.0;        

  // Ad hoc rescaling to try various effects
  adHocSignalrescaler = 1.0;
  adHocDIOrescaler    = 1.0;
  adHocRPCrescaler    = 1.0;
  
  // control of table size and retruned vectors from FofM
  maximumSignalCount = 25;
 
} // setCanonicals


double FMtool::obtainCEdata() 
{
// Returns a number of CE's generated, for use when normalizing.
  // TODO - get the instructions for obtaining this data 
  // from a file or from the user.

  if (OUTPUT_signalEfficiency) {
    std::cout << 
        "\nObtaining conversion electron signal efficiency spectrum...\n";
  }

  // Obtain the CE data as a vector of momenta
  double numberOfGeneratedMuonConversions;
  
  std::vector<std::string> CEfileList;
//#define UseCEFilesMixedForCEdata
#ifdef  UseCEFilesMixedForCEdata
  CEfileList.push_back(
        std::string("/mu2e/data/users/mf/CEFilesMixed.root") );
  numberOfGeneratedMuonConversions = 50000.;
#endif
//#define UseCEFilesMixedAllTime_CE34322ForCEdata
#ifdef  UseCEFilesMixedAllTime_CE34322ForCEdata
  CEfileList.push_back(
     std::string("/mu2e/data/users/mf/CEFilesMixedAllTime/CE34322_0_49.root") );
  numberOfGeneratedMuonConversions = 490000.;
#endif
#define UseCEFilesBrown_CE36825ForCEdata
#ifdef  UseCEFilesBrown_CE36825ForCEdata
  CEfileList.push_back(
     std::string("/mu2e/data/users/mf/CEFilesBrown/CE36825_0_9.root") );
  numberOfGeneratedMuonConversions = 100000.;
#endif
  sigEfficiency.clear();
  sigEfficiency.resize(nBins);
  extractFitmom ( CEfileList, sigEfficiency );
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

  // TODO - provide some for user to provide the # of generated muons

  double POTsRepresented = numberOfGeneratedMuonConversions /
     ( stoppedMuonsPerPOT * capturedMuonsPerStoppedMuon * liveGateFraction );
  double sigEffCountNormalization =  1.0/POTsRepresented;
  double sigEffNormalization =  sigEffCountNormalization / binSize;

  if ( adHocSignalrescaler != 1.0 ) {
    std::cout << "\n**** For this trial, signal normalization is scaled to "
              << adHocSignalrescaler << " times its actual value ****\n\n";
    sigEffNormalization *= adHocSignalrescaler;
  }
 
  if (OUTPUT_signalEfficiency)
  { std::cout << "Signal Efficiency: \n"; }
  unsigned int canonicalRangeCEcount = 0;
  for (unsigned int i=0; i < sigEfficiency.size(); ++i) {
    double p = lowestBin + binSize*i;
    if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) ) {
      canonicalRangeCEcount += sigEfficiency[i];
    }
  }
  if (OUTPUT_signalEfficiency) {
    std::cout << "There are " <<  canonicalRangeCEcount
              << " CE's in the files in the range of "
              << canonicalRangeLo << " - "  << canonicalRangeHi << " \n";
    std::cout << "With normalization, this amounts to "
              << canonicalRangeCEcount * sigEffCountNormalization
                                       * protonsOnTarget * 1.0e-16
               << " CE events in range if BR = 1.0E-16\n";
  }
  for (unsigned int i=0; i < sigEfficiency.size(); ++i) {
    sigEfficiency[i] *= sigEffNormalization;
    if (OUTPUT_signalEfficiency && OUTPUT_spectra) {
      std::cout << lowestBin + binSize*i << ":   " 
                << sigEfficiency[i]/sigEffNormalization;
      std::cout << " --> " << sigEfficiency[i] 
                << " ==> " << sigEfficiency[i] * protonsOnTarget << "\n";
    }
  }
  double integratedSignalEfficiencyInRange = 0;
  for (unsigned int i=0; i < sigEfficiency.size(); ++i) {
    double p = lowestBin + binSize*i ;
    if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) ) {
      integratedSignalEfficiencyInRange += sigEfficiency[i];
    } 
  }
} // normalizeSignalEfficiency

double FMtool::obtainDIOdata(bool & use_diowt) 
{
  // TODO - get the instructions for obtaining this data 
  // from a file or from the user.

  if (OUTPUT_backgroundStrength) {
    std::cout << "\nObtaining DIO spectrum...\n";
  }

  // Obtain the DIO data as a vector of momenta 
  double numberOfGeneratedDIOs;
  
  // Since there will be multiple files for DIO's, we use a list.
  // (Because the TTree is placed under a directory, we cannot use
  // a TChain, we must use a list of TFiles.  Oh well...)
  std::vector<std::string> DIOfileList;
//#define UseDIOFilesMixedntagForDIOdata
#ifdef  UseDIOFilesMixedntagForDIOdata
  DIOfileList.push_back(
              std::string("/mu2e/data/users/mf/DIOFilesMixed1tag.root") );
  DIOfileList.push_back(
              std::string("/mu2e/data/users/mf/DIOFilesMixed2tag.root") );
  numberOfGeneratedDIOs = 194000;
#endif
//#define UseDIOFilesMixedAllTime_DIO34393ForDIOdata
#ifdef  UseDIOFilesMixedAllTime_DIO34393ForDIOdata
  DIOfileList.push_back(
   std::string("/mu2e/data/users/mf/DIOFilesMixedAllTime/DIO34393_0_49.root"));
  DIOfileList.push_back(
   std::string("/mu2e/data/users/mf/DIOFilesMixedAllTime/DIO34393_50_99.root"));
  DIOfileList.push_back(
   std::string("/mu2e/data/users/mf/DIOFilesMixedAllTime/DIO34393_100_199.root"));
  numberOfGeneratedDIOs = 2815000.;
#endif
#define UseDIOFilesBrown_DIO36831ForCEdata
#ifdef  UseDIOFilesBrown_DIO36831ForCEdata
  DIOfileList.push_back(
     std::string("/mu2e/data/users/mf/DIOFilesBrown/DIO36831_0_99.root") );
  DIOfileList.push_back(
     std::string("/mu2e/data/users/mf/DIOFilesBrown/DIO36831_100_199.root") );
  DIOfileList.push_back(
     std::string("/mu2e/data/users/mf/DIOFilesBrown/DIO36831_200_299.root") );
  DIOfileList.push_back(
     std::string("/mu2e/data/users/mf/DIOFilesBrown/DIO36831_300_399.root") );
  numberOfGeneratedDIOs = 10000000.;
#endif

  DIObackground.clear();
  DIObackground.resize(nBins);
  double count103 = 0.0;
  double count104 = 0.0;
  double countRatioThreshold = 5.0;  
  countMcmom ( DIOfileList, count103, count104 ); 
  use_diowt = true;
  if (count104 <= 0) {
    std::cout << "There are no DIO's near 104 in this set.\n"
              << "Apparently, the DIO tracks were generated according to the "
              << "momentum distribution\n" 
              << "and need not be weighted\n";
    use_diowt = false;
  } else if (count103/count104 > countRatioThreshold){
    std::cout << "The ratio of DIO's near 103 to near 104 in this set is "
              << count103/count104 << " -- \n"
              << "Apparently, the DIO tracks were generated according to the "
              << "momentum distribution\n" 
              << "and need not be weighted\n";
    use_diowt = false;
  } else {
    std::cout << "The ratio of DIO's near 103 to near 104 in this set is "
              << count103/count104 << " -- \n"
              << "Apparently, the DIO tracks were generated according to a "
              << "roughly flat momentum distribution.\n" 
              << "The tracks will be weighted using diowt\n";
    use_diowt = true;
  }

  extractFitmom ( DIOfileList, DIObackground, use_diowt );
  
  return numberOfGeneratedDIOs;
  
} // obtainDIOdata

double FMtool::obtainRPCdata() 
{
  // TODO - get the instructions for obtaining this data 
  // from a file or from the user.

  if (OUTPUT_backgroundStrength) {
    std::cout << "\nObtaining RPC spectrum...\n";
  }

  double numberOfGeneratedStoppedPions;
  
  // Obtain the RPC data as a vector of momenta
  // Since there will be multiple files for this, we use a list.
  std::vector<std::string> RPCfileList;

#define UseRPCfiles_RPCFilesMixedxtag_ForRPCdata
#ifdef  UseRPCfiles_RPCFilesMixedxtag_ForRPCdata
  RPCfileList.push_back(
              std::string("/mu2e/data/users/mf/RPCFilesMixed1tag.root") );
  RPCfileList.push_back(
              std::string("/mu2e/data/users/mf/RPCFilesMixed2tag.root") );
  numberOfGeneratedStoppedPions = 2.0e7;
#endif
  RPCbackground.clear();
  RPCbackground.resize(nBins);

  // There is never a menaingful minimum_t0 to cut on in the RPC data 
  // because the distribution drops so rapidly one would never create an honest 
  // one.  Instead, analytic weighting is always used.  So we turn off cutting.
  double remember_minimum_t0 = minimum_t0;
  minimum_t0 = 0;
  extractFitmom ( RPCfileList, RPCbackground );
  minimum_t0 = remember_minimum_t0;

  std::cout << "Enter extinction level, times 10^10: ";
  std::cin >> extinction;
  extinction *= 1.0e-10;
  std::cout << "Using extinction of " << extinction << "\n";

  return numberOfGeneratedStoppedPions;
  
} // obtainRPCdata

void FMtool::extractFitmom 
  ( TTree * tracks, std::vector<double> & counts, bool use_diowt )
{
  // Obtain the momentum data as a vector of momenta
  std::vector<double> momenta;
  std::vector<double> diowts;
  tracks->SetBranchStyle(0);
  int   nactive;   tracks->SetBranchAddress("nactive"  , &nactive);
  float t0err;     tracks->SetBranchAddress("t0err"    , &t0err);    
  float fitmomerr; tracks->SetBranchAddress("fitmomerr", &fitmomerr);    
  float fitcon;    tracks->SetBranchAddress("fitcon"   , &fitcon);    
  int   fitstatus; tracks->SetBranchAddress("fitstatus", &fitstatus);    
  float fitmom;    tracks->SetBranchAddress("fitmom"   , &fitmom);    
  float diowt;     tracks->SetBranchAddress("diowt"    ,  &diowt);    
  float t0;        tracks->SetBranchAddress("t0"       ,  &t0);    
  int ntracks = 0;
  for (int i = 0; true; ++i) {
    int bytesRead = tracks->GetEntry(i);
    if (bytesRead == 0) break;
    ++ntracks;
    if  (    (fitstatus == 1) 
          && (nactive   >= minimum_nactive)
          && (t0err     <= maximum_t0err)
          && (fitmomerr <= maximum_fitmomerr)
          && (fitcon    >= minimum_fitcon) 
          && (fitmom    < 130.0) // deals with the one really screwy track    
        ) 
    {
      if ( minimum_t0 == 0 || t0 > minimum_t0 ) {
        momenta.push_back( fitmom );           
        if ( use_diowt ) {
          diowts.push_back( diowt );
        }
      }
    }     
  } 
  std::cout << "There are " << ntracks << " entries in the tree\n";
  std::cout << "There are " << momenta.size() 
            << " entries passing the cuts\n";
  double binSize = (topOfLastBin - lowestBin)/nBins;
  int nBinnedTracks = 0;
  double dioTotalWeight = 0.0;
  if ( use_diowt ) {
    std::cout << "Extracting DIO's using diowt...\n";
  }
  for (unsigned int i=0; i < momenta.size(); ++i) {
    double p = momenta[i];
    if (p >= lowestBin && p < topOfLastBin) {
      unsigned int binNumber = (p-lowestBin)/binSize;
      if (binNumber >= nBins) binNumber = nBins-1;
      if ( use_diowt ) {
        counts[binNumber] += diowts[i];  // weighted DIO tracks
        dioTotalWeight += diowts[i];
//        if ( (p > 103.5) && (p < 105.0) ) {
//          std::cout << "p = " << p << "  doiwt = " << diowts[i] << "\n";
//        }
      } else {
        counts[binNumber] += 1.0; // everything has equal weight
      }
      ++nBinnedTracks;
    }
  }
  std::cout << "There are " 
            << nBinnedTracks << " tracks between " << lowestBin << " and "
            << topOfLastBin << "\n";
  if ( use_diowt ) {
    std::cout << "(total of DIO weights is " << dioTotalWeight << ")\n";
  }
              
} // extractFitmom



void FMtool::extractFitmom 
  ( std::vector<std::string> const & listOfFileNames
  , std::vector<double> & counts, bool use_diowt )
{
  std::cout << "extracting fitmom from a chain of files \n";
  // Obtain the momentum data as a vector of momenta
  for (unsigned int i=0; i<listOfFileNames.size(); ++i) {
   std::cout << "listOfFileNames["  << i << "].c_str() = "
              << listOfFileNames[i].c_str() << "\n";
    TTreeAccessor ta = accessTTree ( listOfFileNames[i] );
    extractFitmom 
      (  ta.tracks,  counts, use_diowt );
    delete ta.tracks;   ta.tracks = 0;
    delete ta.tfp;      ta.tfp = 0;
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
    std::cout << "(tracks TTree found in ReadKalFits/trkdiag)\n";
  } else {
    ta.tracks = (TTree*)ta.tfp->Get ("trkdiag");
    if (ta.tracks != 0) {
      std::cout << "(tracks TTree found in trkdiag)\n";
    }
  }
  if (ta.tracks == 0) {
    std::cerr << "In file " << fileName << ":\n"
              << "  Failed to open tracks at "  
              << "ReadKalFits/trkdiag or trkdiag" <<  "\n";
    std::exit(1);
  }
  return ta;
} // accessTTree


void FMtool::countMcmom 
  ( TTree * tracks, double & count103, double & count104 )
{
  tracks->SetBranchStyle(0);
  float mcmom;     tracks->SetBranchAddress("mcmom"    , &mcmom);    
  for (int i = 0; true; ++i) {
    int bytesRead = tracks->GetEntry(i);
    if (bytesRead == 0) break;
    if  (  (mcmom >= 102.5) && (mcmom < 103.5) ) {
        count103 += 1.0;    
    } 
    if  (  (mcmom >= 103.5) && (mcmom < 104.5) ) {
        count104 += 1.0;    
    } 
  } 
} // countMcmom

void FMtool::countMcmom 
  ( std::vector<std::string> const & listOfFileNames
  , double & count103, double & count104 )
{
  std::cout << "counting mcmom near 103 and 104 from a chain of files \n";
  count103 = 0.0;
  count104 = 0.0;
  // Obtain the momentum data as a vector of momenta
  for (unsigned int i=0; i<listOfFileNames.size(); ++i) {
    TTreeAccessor ta = accessTTree ( listOfFileNames[i] );
    countMcmom 
      (  ta.tracks,  count103, count104 );
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
  // with p >= 102.5.
  //
  // Sometimes DIO's were generated flat in some range (e.g., 95-105) and
  // should be weighted according to diowt instead.

  // TODO - provide some for user to provide the # of generated DIOs etc.

  double fractionOfDIOsAbove102_5 = 4.42e-15;
  
  double DIOgenerationWeighting = use_diowt ? 1.0 : fractionOfDIOsAbove102_5;

  double DIOsPerStoppedMuon =  1.0 - capturedMuonsPerStoppedMuon;
  double DIObackgroundCountNormalization =
      ( stoppedMuonsPerPOT * DIOsPerStoppedMuon 
          * DIOgenerationWeighting  * liveGateFraction ) /
      ( numberOfGeneratedDIOs );

  if ( adHocDIOrescaler != 1.0 ) {
    std::cout << "\n**** For this trial, DIO normalization is scaled to "
              << adHocDIOrescaler << " times its actual value ****\n\n";
     DIObackgroundCountNormalization *= adHocDIOrescaler;
  }

  DIObackgroundNormalization = DIObackgroundCountNormalization/binSize;
  normalizeAbackground("DIO", DIObackgroundCountNormalization, DIObackground);

} // normalizeDIObackground


void FMtool::normalizeRPCbackground(double numberOfGeneratedRPCs) 
{
  // Normalization of the RPCbackground
  // Note that our convention is to normalize the background to a
  // number per proton on target, so the luminosity does not appear in 
  // the background normalization.

  double internalRPCconversionFactor = 2.0;
  double RPCnormalization =  internalRPCconversionFactor * 9.45e-5;   
       // Explanation:  This is the number to multiply the counts in these 2
        // files by, to get a number of RPC's per 3.6e20 POTs.  It includes
        // the .021 RPC's per stopped pion, and the live gate fraction,
        // which is very small (at t = 792 it is 2.6e-16) because most stopped 
        // pions are in the first 260 nsec, and the lifetime is only 21 nsec.  
        // Since the files in mixed1 and mixed2 were based on generating 2*10^7 
        // stopped pions, this product gives 9.45e-5.  
        //
        // Since we want a number per POT, we'll need to divide that back out.  
        //
        // The number 9.45e-5 is based on just the processes Bob B. had
        // simulated to get these files.  But in fact there are other processes,
        // in particular, internal conversion is almost exactly as intense
        // as as RPC photon converts in the stopping target (though the angular
        // distribution will be a bit different.    
        // An email (rhbob) has said to just multiply the normalization by 2, 
        // (nternalRPCconversionFactor) which we do above.

  if ( minimum_t0 != 0.0 ) {
    double RPCliveGateFraction;
    if (minimum_t0 < 400) {
      RPCliveGateFraction = RPCtimeProbability(400);
    } else if (minimum_t0 > 1000) {
      RPCliveGateFraction = 0;
    } else {
      RPCliveGateFraction = RPCtimeProbability(minimum_t0);
    }
    RPCliveGateFraction += extinction * stoppedPionsPerPOT;
    
    RPCnormalization = protonsOnTarget * RPCperStoppedPion * RPCliveGateFraction
                     * internalRPCconversionFactor / numberOfGeneratedRPCs;  
    std::cout << "with minimum_t0 = " << minimum_t0 
              << " and extinction = " << extinction << "\n"
              << "RPCliveGateFraction is "  
              << RPCliveGateFraction - extinction * stoppedPionsPerPOT
              << " + " << extinction * stoppedPionsPerPOT << "\n"
              << "RPCnormalization is " << RPCnormalization << "\n";
  }
  
  
  double RPCbackgroundCountNormalization =
        ( RPCnormalization ) /
        ( protonsOnTarget );

  if ( adHocRPCrescaler != 1.0 ) {
    std::cout << "\n**** For this trial, RPC normalization is scaled to "
              << adHocRPCrescaler << " times its actual value ****\n\n";
     RPCbackgroundCountNormalization *= adHocRPCrescaler;
  }

  RPCbackgroundNormalization = RPCbackgroundCountNormalization/binSize;
  normalizeAbackground("RPC", RPCbackgroundCountNormalization, RPCbackground);

} // normalizeRPCbackground

void FMtool::normalizeAbackground(  
        std::string const & bname, 
        double backgroundCountNormalization, 
        std::vector<double> & background)
{        
  double backgroundNormalization =  backgroundCountNormalization / binSize;  
             
  if (OUTPUT_backgroundStrength) {
    std::cout << "\n" << bname << "background: (Count normalization is " 
              << backgroundCountNormalization << " ==> " 
              << backgroundCountNormalization*protonsOnTarget << ")\n\n";
  }
  double canonicalRangeCount = 0;
  for (unsigned int i=0; i < background.size(); ++i) {
    double p = lowestBin + binSize*i;
    if ( p >= 103.51 && p <= 104.70 )  {
      canonicalRangeCount += background[i];
    }
  }
  if (OUTPUT_backgroundStrength) {
    std::cout << "There are " <<  canonicalRangeCount
              << " " << bname << "'s in the files in the range of "
              << canonicalRangeLo << " - "  << canonicalRangeHi << " \n";
    std::cout << "With normalization, this amounts to "
              << canonicalRangeCount * backgroundCountNormalization 
                                        * protonsOnTarget 
               << " " << bname << " events in the range \n";
  }
  for (unsigned int i=0; i < background.size(); ++i) {
    background[i] *= backgroundNormalization;
    if (OUTPUT_backgroundStrength && OUTPUT_spectra) {
      std::cout << lowestBin + binSize*i << ":   "
                << background[i] << " --> " 
                << background[i] / backgroundCountNormalization;
      std::cout << " --> " << background[i] 
                << " ==> " << background[i] * protonsOnTarget << "\n";
    }
  }

} // normalizeAbackground


void FMtool::applyFofM() const
{
  //
  // Use the FofM class to calculate quantities of interest
  //

  // Since our vectors represent values at the middle of intervals,
  // we should state our sample points vectors as such when using
  // them to construct spectra
  double midBinOffset = (topOfLastBin - lowestBin)/(2.0*nBins);
  double lowestPoint = lowestBin + midBinOffset;
  double endingPoint = topOfLastBin - midBinOffset;

  std::cout << "\n------------------------------------ \n\n"
            << "Using the FofM tool \n\n";


  std::cout << "lowestBin = " << lowestBin
            << "  topOfLastBin = " << topOfLastBin 
            << "\nlowestPoint = " << lowestPoint
            << "  endingPoint = " << endingPoint
            << "  nBins = " << nBins << "\n\n";    

  // TODO - think out and implement the way we want users to exercise this
  //        control
  double lowCut;
  double highCut;
  std::cout <<
    "Enter fixed low momentum cut (or 0 for automatic optimizing): ";
  std::cin >> lowCut;
  std::cout <<
    "Enter fixed high momentum cut (or 0 for automatic optimizing): ";
  std::cin >> highCut;
  MeritFunctionChoice mfc = SmoothedPunziMeritFunction; 
  int mfcN;
  if ( (lowCut == 0) || (highCut == 0) ) {
    std::cout << "Momentum cut optimization: \n"
              << "Enter 0 to use smoothed-worst-case merit function \n" 
              << "Enter 1 to use 90% CL sensitivity merit function:   ";
    std::cin>> mfcN;
    if (mfcN == 1) mfc = FCsensitivityMeritFunction;
  } 
  FofM::Spectrum CE_spectrum (sigEfficiency, lowestPoint, endingPoint);
  FofM::Spectrum DIO_spectrum (DIObackground, lowestPoint, endingPoint);
  FofM figureOfMeritCalculator ( CE_spectrum, 
                                 DIO_spectrum, 
                                 protonsOnTarget,
                                 mfc );  

  if (OUTPUT_backgroundSplines) { 
    figureOfMeritCalculator.displayBackground(lowestBin, topOfLastBin, 100);
    splines::Spline<1> sbkg = figureOfMeritCalculator.getBackground();
    splines::Grid<1> sgrid = figureOfMeritCalculator.getGrid();
    std::cout << " L * DIO Background (normalized) vs L * Background \n";
    double roughInt = 0.0;
    double gridstep = sgrid[1] - sgrid[0];
    for (unsigned int i = 0; i < sgrid.nPoints(); ++i) {
      double p =  sgrid[i];
      std::cout << i << ": p = " << p  << "  " 
                << DIObackground[i] / DIObackgroundNormalization << "  "
                << protonsOnTarget * DIObackground[i] << "  "
                << protonsOnTarget * sbkg(p) << "\n";
      if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) )  {
        roughInt += sbkg(p);
      }
    } 
    std::cout << "Rough integral in canonical range based on sbkg = " 
              << protonsOnTarget*gridstep*roughInt << "\n";
  }

  FofM::Spectrum RPC_spectrum (RPCbackground, lowestPoint, endingPoint);
  figureOfMeritCalculator.addBackground (RPC_spectrum);

  if (OUTPUT_backgroundSplines) { 
    figureOfMeritCalculator.displayBackground(lowestBin, topOfLastBin, 100);
    splines::Spline<1> sbkg = figureOfMeritCalculator.getBackground();
    splines::Grid<1> sgrid = figureOfMeritCalculator.getGrid();
    std::cout << " L * (DIO+RPC) Background (normalized) vs L * Background \n";
    double roughInt = 0.0;
    double gridstep = sgrid[1] - sgrid[0];
    for (unsigned int i = 0; i < sgrid.nPoints(); ++i) {
      double p =  sgrid[i];
      std::cout << i << ": p = " << p  << "  " 
                << DIObackground[i] / DIObackgroundNormalization 
                +  RPCbackground[i] / RPCbackgroundNormalization << "  "
                << protonsOnTarget * 
                        (DIObackground[i] + RPCbackground[i]) << "  "
                << protonsOnTarget * sbkg(p) << "\n";
      if ( (p >= canonicalRangeLo) && (p <= canonicalRangeHi) )  {
        roughInt += sbkg(p);
      }
    } 
    std::cout << "DIO + RPC Rough integral in canonical range based on sbkg = " 
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

  std::string t; 
  if ((lowCut == 0) && (highCut == 0)) {    
    t = figureOfMeritCalculator.tables(maximumSignalCount);
  } else if (lowCut == 0) { 
    t = figureOfMeritCalculator.tables_fixed_highCut 
                              (maximumSignalCount, highCut);
  } else if (highCut == 0) { 
    t = figureOfMeritCalculator.tables_fixed_lowCut 
                              (maximumSignalCount, lowCut);
  } else {
    t = figureOfMeritCalculator.tables_fixed_cuts 
                              (maximumSignalCount, lowCut, highCut);
  }

  std::cout << t << "\n";

  std::cout << "\n--------------------------------------- \n\n";

} // applyFofM()


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

} // end namespace mu2e

int main () {
  mu2e::FMtool f;
  f.analyze();
  return 0;
}
