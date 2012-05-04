#ifndef FMTOOL_H
#define FMTOOL_H

// ----------------------------------------------------------------------
//
// FMtool_H.h
//
// Tool to use FofM class for calculations based on data from simulations.
// sense of the Cousins/Feldman paper 
//
// ----------------------------------------------------------------------

#include <vector>
#include <string>

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TTree.h"

// mu2e includes
#include "splines/Spline.h"

namespace mu2e {

class FMtool
{
public:
  FMtool();
  void   analyze();
  void   setFilterLevels();
  void   setCanonicals();
  double obtainCEdata();
  void   normalizeSignalEfficiency(double numberOfGeneratedMuonConversions);
  double obtainDIOdata(bool & use_diowt);
  void   normalizeDIObackground(double numberOfGeneratedDIOs, bool use_diowt);  
  double obtainRPCdata();
  void   normalizeRPCbackground(double numberOfGeneratedRPCs);  
  void   normalizeAbackground(  
          std::string const & bname, 
          double backgroundCountNormalization, 
          std::vector<double> & background);
  void   applyFofM() const;
  
private:

  enum TrackerCutsSet {
      TrackerCuts_A
    , TrackerCuts_B
    , TrackerCuts_C
    , TrackerCuts_D 
  };

  void ourTrackerCuts ( TrackerCutsSet cuts, 
                        int & nactive, 
                        float & t0err, 
                        float & fitmomerr,
                        float & fitcon );
  void nextTrackerCuts( TrackerCutsSet & cuts );  
  void extractFitmom 
        ( TTree * tracks
        , std::vector<double> & counts
        , bool use_diowt = false );
  void extractFitmom 
       ( std::vector<std::string> const & listOfFileNames
       , std::vector<double> & counts
       , bool use_diowt = false  );
  struct TTreeAccessor {
    TFile* tfp;
    TTree* tracks;
  };
  TTreeAccessor accessTTree (std::string const & fileName) const;
  void countMcmom ( TTree * tracks, double & count103, double & count104);
  void countMcmom ( std::vector<std::string> const & listOfFileNames, 
                        double & count103, double & count104);

  std::vector<double> RPCtimeProbabilityVector() const;
  splines::Spline<1> RPCtimeProbabilitySpline()  const;

public:
  bool OUTPUT_signalEfficiency;
  bool OUTPUT_spectra;
  bool OUTPUT_backgroundStrength;
  bool OUTPUT_backgroundSplines;

  double adHocSignalrescaler;
  double adHocDIOrescaler;
  double adHocRPCrescaler;

private:

  // filter information
  int    minimum_nactive;  
  float  maximum_t0err;    
  float  maximum_fitmomerr;
  float  minimum_fitcon;   
  float  minimum_t0;
    
  // canonical numbers
  double canonicalRangeLo;
  double canonicalRangeHi;

  // table control
  int    maximumSignalCount;
  
  // binning information
  double lowestBin;  
  double topOfLastBin; 
  unsigned int nBins;
  double binSize;

  // Experiment quantities
  double protonsOnTarget;  
  double stoppedMuonsPerPOT;
  double capturedMuonsPerStoppedMuon;
  double liveGateFraction;
  double RPCperStoppedPion;
    
  // Spectra
  std::vector<double> sigEfficiency;
  std::vector<double> DIObackground;
  std::vector<double> RPCbackground;
  
  // normalization factors used
  double DIObackgroundNormalization;
  double RPCbackgroundNormalization;

  // time window shapes
  splines::Spline<1> RPCtimeProbability;
  double stoppedPionsPerPOT;
  double extinction;
   
}; // FMtool

} // end namespace mu2e

#endif // FMTOOL_H
