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

  ParameterSet cuts = pset.get<ParameterSet>("trackerCuts");

  minimum_nactive   = cuts.get<int>("minimum_nactive"); 
  maximum_t0err     = cuts.get<float>("maximum_t0err"); 
  maximum_fitmomerr = cuts.get<float>("maximum_fitmomerr"); 
  minimum_fitcon    = cuts.get<float>("minimum_fitcon"); 
  
  tCuts = pset.get< std::vector<double> >("t0cuts");

}  // setFilterLevels


void FMtool::setCanonicals() 
{
  // TODO - allow this to be set from file or whatever.

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
  cadence = 1694;  // nsec after start of proton pulse
  stoppedMuonsPerPOT = 0.0021;  // I have seen .00215; new  number is .0016
  capturedMuonsPerStoppedMuon = 0.609; // DocDB 48 - I have seen .59
  RPCperStoppedPion = 0.021;     // DocDB 1087
  stoppedPionsPerPOT = 1.53e-6;
  
  // control of table size and returned vectors from FofM
  maximumSignalCount = 25;
 
} // setCanonicals


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
     
  size_t nTcuts = tCuts.size();
  sigEfficiency.clear();
  sigEfficiency.resize(nTcuts);
  for (size_t tCutNumber = 0; tCutNumber < nTcuts; ++tCutNumber) {
    sigEfficiency[tCutNumber].clear();
    sigEfficiency[tCutNumber].resize(nBins);
  }
  extractFitmom ( CEfileList, sigEfficiency );

  if (OUTPUT_dataProperties) {
    os << "CE's based on " << numberOfGeneratedMuonConversions 
       << " generated muon conversions \n"; 
    if (CEliveGateFraction != 1.0) {
      os << "CE data used a live time gate fraction of " 
         << CEliveGateFraction << "\n";
    }
  }
  adHocSignalrescaler = CEpset.get<double>("adHocSignalrescaler",1.0);

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

  os << "sigEffNormalization = " << sigEffNormalization << "\n";

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

  double numberOfGeneratedDIOs = DIOpset.get<double>("numberOfGeneratedDIOs");

  DIOliveGateFraction = DIOpset.get<double>("liveGateFraction");
    
  DIObackground.clear();
  DIObackground.resize(nBins);
  double count103 = 0.0;
  double count104 = 0.0;
  double countRatioThreshold = 5.0;  
  countMcmom ( DIOfileList, count103, count104 ); 

  use_diowt = DIOpset.get<bool>("dioTracksWeighted");

  os << "The ratio of DIO's near 103 to near 104 in this set is "
     << count103/count104 << " -- \n";

  if (count104 <= 0) {
    if ( use_diowt ) {
      std::cerr << "Warning: Discrepancy detected involving use of diowt! \n"; 
      os << "There are no DIO's near 104 in this set.\n"
         << "Apparently, the DIO tracks were generated \n"
         << "according to the DIO momentum distribution " 
         << "and should not not be weighted.\n"
         << "But dioTracksWeighted (in the fcl) has value true.\n";
      os        << "Warning: Discrepancy detected involving use of diowt! \n"; 
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
    }
  }
 
  if ( use_diowt ) {
    fractionOfDIOsRepresented = 1.0;
  } else {
    fractionOfDIOsRepresented = 
      DIOpset.get<double>("fractionOfDIOsRepresented");
  }

  size_t nTcuts = tCuts.size();
  DIObackground.clear();
  DIObackground.resize(nTcuts);
  for (size_t tCutNumber = 0; tCutNumber < nTcuts; ++tCutNumber) {
    DIObackground[tCutNumber].clear();
    DIObackground[tCutNumber].resize(nBins);
  }
  extractFitmom ( DIOfileList, DIObackground, use_diowt );
  
  if (OUTPUT_dataProperties) {
    os << "DIO's based on " << numberOfGeneratedDIOs 
       << "  generated DIO electrons \n"; 
    if (DIOliveGateFraction != 1.0) {
      os << "DIO data used a live time gate fraction of " 
         << DIOliveGateFraction << "\n";
    }
  }
  adHocDIOrescaler = DIOpset.get<double>("adHocDIOrescaler",1.0);

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
  
  return numberOfGeneratedRPCs;
  
} // obtainRPCdata

void FMtool::extractFitmom 
  ( TTree * tracks, 
    std::vector< std::vector<double> > & counts, 
    bool use_diowt, 
    bool apply_timeCuts )
{
  // Set up place to obtain the momentum and other data for immediate use 
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

  double binSize = (topOfLastBin - lowestBin)/nBins;
  std::vector<int> nBinnedTracks(tCuts.size());
  std::vector<double> dioTotalWeight(tCuts.size());
  size_t nt = tCuts.size();
  for (int i = 0; true; ++i) {
    int bytesRead = tracks->GetEntry(i);
    if (bytesRead == 0) break;
    ++ntracks;
    if  (    (fitstatus == 1)     // track quality cuts 
          && (nactive   >= minimum_nactive)
          && (t0err     <= maximum_t0err)
          && (fitmomerr <= maximum_fitmomerr)
          && (fitcon    >= minimum_fitcon) 
          && (fitmom    < 120.0) // deals with any really screwy tracks    
        ) 
    {
      for (size_t t = 0; t < nt; ++t) {
        if ( (!apply_timeCuts) || (t0 >= tCuts[t]) ) {   // found a good track!
          double p = fitmom;
          if (p >= lowestBin && p < topOfLastBin) { // p in range to be binned
            unsigned int binNumber = (p-lowestBin)/binSize;
            if (binNumber >= nBins) binNumber = nBins-1;
            if ( use_diowt ) {
              counts[t][binNumber] += diowt;  // weighted DIO tracks
              dioTotalWeight[t] += diowt;
            } else {
              counts[t][binNumber] += 1.0; // everything has equal weight
            }
            ++nBinnedTracks[t];
          } // end of *if* on p in range to be binned
        } // end of *if* on found a good track
      } // end of *for* on t going up to nt 
    } // end of *if* checking the track quality cuts    
  } // end of *for* loop on i which goes thru track entries
  if (OUTPUT_fileEntriesStatistics) {
    os << "There are " << ntracks << " entries in the tree\n";
    for ( size_t i = 0; i < nt; ++i ) {
      os << "There are " << nBinnedTracks[i] 
         << " binned entries passing the cuts with t0 >= " 
         << tCuts[i] << "\n";
      if (use_diowt) { 
        os << "(total of DIO weights passing cuts with t0 >= "
           <<  tCuts[i] << " is " << dioTotalWeight[i] << ")\n";
      }
    }
  }
              
} // extractFitmom

void FMtool::extractFitmom 
  ( std::vector<std::string> const & listOfFileNames
  , std::vector< std::vector<double> > & counts
  , bool use_diowt
  , bool apply_timeCuts 
  )
{
  // Obtain the momentum data as a vector of momenta
  if (OUTPUT_fileNames) {
    os << "extracting fitmom from a chain of files \n";
  }
  for (unsigned int i=0; i<listOfFileNames.size(); ++i) {
    if (OUTPUT_fileNames) {
      os << "["  << i << "] : "
                << listOfFileNames[i].c_str() << "\n";
    }
    TTreeAccessor ta = accessTTree ( listOfFileNames[i] );
    extractFitmom 
      (  ta.tracks,  counts, use_diowt, apply_timeCuts );
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
  // with p >= 102.5 (or, in principle, some other lower limit).  If
  // the tracks were generated with a momentum distribution following
  // czarnecki, then the tracks should not be weghted, but an overall 
  // factor representing the fraction of DIO's above the minimum momnetum
  // must be used.  For a minimum p of 102.5, this factor is 4.42e-15.
  //
  // In other files, DIO's were generated flat in some momentum range 
  // (e.g., 95-105) in which case the tracks should be weighted according 
  // to diowt.  The typical diowt is in the area of 1.0e-15.

  double DIOgenerationWeighting = use_diowt ? 1.0 : fractionOfDIOsRepresented;

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

  DIObackgroundNormalization = DIObackgroundCountNormalization/binSize;
  for (size_t t =0; t <  tCuts.size(); ++t) {
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
  // We start by forming RPCnormalization which is the number of RMC
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

  ParameterSet RPCpset = 
    pset.get<ParameterSet>("radiativePionCaptureElectrons");  
  
  double RPCnormalizationBase =  RPCperStoppedPion / numberOfGeneratedRPCs;  
  // RPCnormalizationBase still needs a factor of RPCtimeProbability

  if ( adHocRPCrescaler != 1.0 ) {
     os << "\n**** For this trial, RPC normalization is scaled to "
              << adHocRPCrescaler << " times its actual value ****\n\n";
     RPCnormalizationBase *= adHocRPCrescaler;
  }

  for (size_t t = 0; t <  tCuts.size(); ++t) {
    double RPClivePionFraction =  RPCtimeProbability(tCuts[t]);
    double flatPionBackgroundInLiveGate = extinction * stoppedPionsPerPOT 
                * (cadence - tCuts[t]) / cadence;
        // This is the number of stopped pions arriving in the live gate
        // (between the cut time and the end of the tiem window) per
        // desirable POT, due to the (tiny) number of POT's outside the
        // proton pulse (the non-extinct protons)   
    if (OUTPUT_RPClivePionFraction) {
      os << "RPClivePionFraction = " << RPClivePionFraction << " + " 
         << flatPionBackgroundInLiveGate << " = " 
         << RPClivePionFraction + flatPionBackgroundInLiveGate << "\n";
    }
    
    RPClivePionFraction += flatPionBackgroundInLiveGate;
    double RPCbackgroundCountNormalization = 
                RPCnormalizationBase * RPClivePionFraction;
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
    if ( p >= 103.51 && p <= 104.70 )  {
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

void FMtool::applyFofM() const
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

  os << "\n------------------------------------ \n\n"
     << "Using the FofM tool \n\n";


  os << "lowestBin = " << lowestBin
     << "  topOfLastBin = " << topOfLastBin 
     << "\nlowestPoint = " << lowestPoint
     << "  endingPoint = " << endingPoint
     << "  nBins = " << nBins << "\n\n";    

  // TODO - think out and implement the way we want users to exercise this
  //        control
  double lowCut;
  double highCut;
  ParameterSet pCutpset = pset.get<ParameterSet>("momentumCuts");
  bool optimize_pcut_low = pCutpset.get<bool>("optimize_pcut_low",true);
  lowCut = optimize_pcut_low ? 0 : pCutpset.get<double>("fixed_pcut_low");
  bool optimize_pcut_high = pCutpset.get<bool>("optimize_pcut_high",true);
  highCut = optimize_pcut_high ? 0 : pCutpset.get<double>("fixed_pcut_high");
  
  bool useSmoothedPunziFMtoOptimize = pset.get<bool>
                        ("useSmoothedPunziFMtoOptimize", false);
                         
  MeritFunctionChoice mfc = useSmoothedPunziFMtoOptimize ? 
        SmoothedPunziMeritFunction : FCsensitivityMeritFunction; 

  std::string table;
  std::vector<FofM::Summary> summaries(tCuts.size());
  double best_tCutMerit  = 0;
  size_t best_tCutNumber = 0;
  for (size_t tCutNumber=0; tCutNumber<tCuts.size(); ++tCutNumber) {
    summaries[tCutNumber] = 
        applyFofM( tCutNumber, lowestPoint, endingPoint, table, 
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
     << " pLow     Phigh \n";
  for (size_t tCutNumber=0; tCutNumber<tCuts.size(); ++tCutNumber) {
    os << tCuts[tCutNumber] 
       << "  " << summaries[tCutNumber].figureOfMerit;
    if (tCutNumber == best_tCutNumber) {
      os << " ** ";
    } else {
      os << "    ";
    }
    os <<  summaries[tCutNumber].singleEventSensitivity 
       << "  " << summaries[tCutNumber].CL90sensitivity 
       << "  " << summaries[tCutNumber].smoothedPunziSensitivity 
       << "  "  <<  summaries[tCutNumber].pCutLo
       << "  "  <<  summaries[tCutNumber].pCutHi << "\n";
  }

// Repeat the best one if there were more than 2 t cuts explored:
  if ( (tCuts.size() > 2) || (!OUTPUT_allTables) ) {
    summaries[best_tCutNumber] = 
        applyFofM( best_tCutNumber, lowestPoint, endingPoint, table, 
                        mfc, lowCut, highCut );
    os << table << "\n--------------------------------------- \n\n";
  }
  
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

  FofM::Spectrum CE_spectrum (sigEfficiency[tCutNumber], 
                                lowestPoint, endingPoint);
  FofM::Spectrum DIO_spectrum (DIObackground[tCutNumber], 
                                lowestPoint, endingPoint);
  FofM figureOfMeritCalculator ( CE_spectrum, 
                                 DIO_spectrum, 
                                 protonsOnTarget,
                                 mfc );  

  if (OUTPUT_backgroundSplines) { 
    figureOfMeritCalculator.displayBackground(lowestBin, topOfLastBin, 100);
    splines::Spline<1> sbkg = figureOfMeritCalculator.getBackground();
    splines::Grid<1> sgrid = figureOfMeritCalculator.getGrid();
    os << " L * DIO Background (normalized) vs L * Background \n";
    double roughInt = 0.0;
    double gridstep = sgrid[1] - sgrid[0];
    for (unsigned int i = 0; i < sgrid.nPoints(); ++i) {
      double p =  sgrid[i];
      os << i << ": p = " << p  << "  " 
         << DIObackground[tCutNumber][i] / DIObackgroundNormalization << "  "
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
  FofM::Spectrum RPC_spectrum (RPCbackground[tCutNumber], 
                                lowestPoint, endingPoint);
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
    t = figureOfMeritCalculator.tables(maximumSignalCount, os, summary);
  } else if (lowCut == 0) { 
    t = figureOfMeritCalculator.tables_fixed_highCut 
                              (maximumSignalCount, highCut, os, summary);
  } else if (highCut == 0) { 
    t = figureOfMeritCalculator.tables_fixed_lowCut 
                              (maximumSignalCount, lowCut, os, summary);
  } else {
    t = figureOfMeritCalculator.tables_fixed_cuts 
                            (maximumSignalCount, lowCut, highCut, os, summary);
  }

  os << t << "\n";

  os << "\n--------------------------------------- \n\n";

  return summary;
  
} // applyFofM(cutNumber...)


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

  std::string outputName = pset.get<std::string>("outputFile"); 
  if (outputName == "cout") {
    mu2e::FMtool f (pset, std::cout) ;
    f.analyze();
  } else {
    std::ofstream os (outputName.c_str());
    mu2e::FMtool f (pset, os) ;
    f.analyze();
  }
  return 0;
}
