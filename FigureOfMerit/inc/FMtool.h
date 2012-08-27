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
#include <iostream>

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TTree.h"

// mu2e includes
#include "splines/Spline.h"
#include "FigureOfMerit/inc/FofM.h"

// art and art externals includes
#include "fhiclcpp/ParameterSet.h"

namespace mu2e {

class FMtool
{
public:
  FMtool (fhicl::ParameterSet const & p, std::ostream & o);
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
  FofM::Summary applyFofM( 
                    size_t tCutNumber,  double lowestPoint,
                    double endingPoint, std::string & table,
                    MeritFunctionChoice mfc,
                    double lowCut, double highCut ) const;
  
private:
  class NthHighest;
  
  void decideVerbosity();
  void extractFitmom 
        ( TTree * tracks
        , std::vector< std::vector<double> > & counts
        , std::vector<NthHighest> & lowPs
        , bool dio = false
        , bool use_diowt = false, bool apply_timeCuts = true );
  void extractFitmom 
       ( std::vector<std::string> const & listOfFileNames
       , std::vector< std::vector<double> > & counts, bool dio = false
       , bool use_diowt = false, bool apply_timeCuts = true  );
  struct TTreeAccessor {
    TFile* tfp;
    TTree* tracks;
  };
  TTreeAccessor accessTTree (std::string const & fileName) const;
  void countMcmom ( TTree * tracks, 
                        double & count103, double & count104, bool use_diowt);
  void countMcmom ( std::vector<std::string> const & listOfFileNames, 
                        double & count103, double & count104, bool use_diowt);

  std::vector<double> RPCtimeProbabilityVector() const;
  splines::Spline<1> RPCtimeProbabilitySpline()  const;

  void DIOstatisticsWarning
                (double computedLowCut, size_t tCutNumber) const;
  void DIOstatisticsWarningFixedLowCut
                (double fixedLowCut, size_t tCutNumber) const;
public:
  // hard-coded verbosity flags

  // fcl controlled verbocity flags
  bool OUTPUT_signalEfficiency;
  bool OUTPUT_spectra;
  bool OUTPUT_backgroundStrength;
  bool OUTPUT_backgroundSplines;
  bool OUTPUT_progress;
  bool OUTPUT_tracksTTreeLocation;
  bool OUTPUT_fileNames;
  bool OUTPUT_fileEntriesStatistics;
  bool OUTPUT_dataProperties;
  bool OUTPUT_RPClivePionFraction;
  bool OUTPUT_allTables;
  bool OUTPUT_detailedTrace;
  bool OUTPUT_detailedCutTrace;
  
  double adHocSignalrescaler;
  double adHocDIOrescaler;
  double adHocRPCrescaler;

private:

  // parameter input
  fhicl::ParameterSet const & pset;

  // output file
  std::ostream & os;

  // filter information
  int    minimum_nactive;  
  float  maximum_t0err;    
  float  maximum_fitmomerr;
  float  minimum_fitcon;   
  float  minimum_t0;

  // time windows
  std::vector<double> tCuts;
  double  CEliveGateFraction;
  double DIOliveGateFraction;   
    
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
  double cadence;  // inverse frequency, nsec
  double liveGateFraction;
  double RPCperStoppedPion;
    
  // Spectra
  std::vector< std::vector<double> > sigEfficiency;
  std::vector< std::vector<double> > DIObackground;
  std::vector< std::vector<double> > RPCbackground;
  
  // normalization factors used
  double fractionOfDIOsRepresented;
  double DIObackgroundNormalization;
  double RPCbackgroundNormalization;
  double DIOflatGenerationWindowLo;
  double DIOflatGenerationWindowHi;
  
  // input counting
  int tracksExtracted;
  int tracksBinned;
  
  // Protection agains low-statistics DIO tails
  int minimumMeaningfulDIOtail;
  std::vector<double> lowMomentumCutCeiling;
  class NthHighest {
  public:
    NthHighest (int nn, double topp ) 
    : n(nn)
    , m(0)
    , top(topp)
    , vals(n,top) {}
    NthHighest () : n(0) , m(0) {}
    void add (double p);
    double value () const {return (m>0 ? vals[0] : 0.0);}
  private:
    int n;
    int m;
    double top;
    std::vector<double> vals;
  };

    
  // time window shapes
  splines::Spline<1> RPCtimeProbability;
  double stoppedPionsPerPOT;
  double extinction;
  
  
}; // FMtool

} // end namespace mu2e

#endif // FMTOOL_H
