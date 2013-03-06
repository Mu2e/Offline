#ifndef FMTOOL_H
#define FMTOOL_H

// ----------------------------------------------------------------------
//
// FMtool_H.h
//
// Tool to use FofM class for calculations based on data from simulations.
//
// ----------------------------------------------------------------------

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

// Root includes.
#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TVectorT.h"

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
  void   setRootGraphics();
  void   establishHistograms();
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
  void   applyFofM();
  FofM::Summary applyFofM( 
                    size_t tCutNumber,  double lowestPoint,
                    double endingPoint, std::string & table,
                    MeritFunctionChoice mfc,
                    double lowCut, double highCut ) const;
  
private:
  class NthHighest;
  class RootGraphics;

  enum {
        CEfitmom
      , DIOfitmom
      , RPCfitmom
  } SpectrumType;
    
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
  void accumulateCutStatistics( int fitstatus, int nactive,  double fitcon,
				double fitmomerr, double t0err, double fitmom, 
			        int & nfitstatus,  int & nnactive, 
			        int & nfitcon, int & nfitmomerr, 
				int & nt0err, int & nfitmom ) const; 
  void outputCutSequence 
    ( int treeEntries,int nfitstatus,int nnactive,
      int nfitcon, int nfitmomerr, int nt0err, int nfitmom ) const; 
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
  void DIOnegativeExplanation 
                (double fixedLowCut, size_t tCutNumber, double newFixedLowCut) const;

  static std::vector<double> rebin(std::vector<double> const & x, unsigned int n, 
				   double a, double b);
  void closeRootFiles();

  void recordFclGroup(std::string const & name) const {
    os << "                    " << name << "\n";
  }
  template<typename T> void recordFcl(std::string const & name, T const & val) const {
    os << "                        " << name << " : " << val << "\n";
  }

  void bottomLine ( std::string const & description, double value );

  void makeDIOreferenceTF1(int tracksExtracted, double DIOflatGenerationWindowLo);

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
  unsigned int nSplineBins;
  
  // Experiment quantities
  double protonsOnTarget;  
  double stoppedMuonsPerPOTold;
  double stoppedMuonsPerPOT;
  double capturedMuonsPerStoppedMuon;
  double cadence;  // inverse frequency, nsec
  double liveGateFraction;
  double RPCperStoppedPion;
    
  int stoppedMuonsDef;
  int stoppedMuonsThisRun;
  double sc_factor;


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
    double value () const {return (m>0 ? vals[0] : top);}
    int howManyHigher(double p) const;
  private:
    int n;
    int m;
    double top;
    std::vector<double> vals;
  };
  std::vector<NthHighest> highPs;
    
  // time window shapes
  splines::Spline<1> RPCtimeProbability;
  double stoppedPionsPerPOT;
  double extinction;

  // "bottom line" information
  std::ostringstream botMenu;
  std::ostringstream botLine;  

  //public:
  //float val0; 

  // Root graphics 
  struct RootGraphics
  {
    bool enabled;
    std::string rootFile;
    std::string cintFile;
    TFile* rootfp;
    TTree* tree;
    TBranch* branch;
    std::ostringstream cint; 

    static const int nvalv = 20;    float valv[nvalv];
    static const int nmomv = 5;     float momv[nmomv];

    static const int nhist = 3;  
    TH1F *h_fitmomCE[nhist];
    TH1F *h_fitmomDIO[nhist]; 
    TH1F *hQ_fitmomCE[nhist];
    TH1F *hQ_fitmomDIO[nhist];


    TH1F* hCEspectrum;
    TH1F* hCEspectrumX;
    TH1F* hDIOspectrum;
    TH1F* hDIOspectrumX;

    TF1*  fDIOreference;





    RootGraphics() 
      : enabled(false)
      , rootfp(0)
      , tree(0)
      , branch(0)
      , hDIOspectrum(0)
      , hDIOspectrumX(0)
      , fDIOreference(0)
    {}
    ~RootGraphics() { if (enabled) delete rootfp; }
    bool write_cint_file() const;   
  }; 

  RootGraphics rg;
}; // FMtool


} // end namespace mu2e

#endif // FMTOOL_H
