///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CalPatRec_CalTimePeakFinder_module
#define CalPatRec_CalTimePeakFinder_module

#ifdef __GCCXML__A
namespace art {
  //  class EDProducer;
  class EDFilter;
  class Run;
  class Event;
};
#else
#  include "art/Framework/Core/EDFilter.h"
#  include "art/Framework/Principal/Event.h"
#endif

// data
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/HelixVal.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/TimeClusterCollection.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BField/BField.hh"
#include "BTrk/BField/BFieldFixed.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "Mu2eBTrk/inc/BaBarMu2eField.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "BTrk/BaBar/BbrStringUtils.hh"
#include "CalPatRec/inc/TrkDefHack.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"
#include "CalPatRec/inc/CalTimePeak.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"

#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/Doublet.hh"

#include "TROOT.h"
#include "TFolder.h"
#include "CalPatRec/inc/KalFitHack.hh"
#include "CalPatRec/inc/HelixFitHack.hh"
#include "CalPatRec/inc/THackData.hh"

// Mu2e
#include "RecoDataProducts/inc/KalRepPayloadCollection.hh"
#include "TrkPatRec/inc/PayloadSaver.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"

//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TGMsgBox.h"
#include "TTree.h"
#include "TFolder.h"

//#include "TStopwatch.h"
// #include "TSpectrum.h"
// #include "TSpectrum2.h"
// #include "TSpectrum3.h"
// #include "TMVA/Reader.h"
// boost
// C++
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <functional>
#include <float.h>
#include <vector>
#include <set>
#include <map>

class Ref;
class THackData;

namespace fhicl {
  class ParameterSet;
}



namespace mu2e {  
  class Calorimeter;
  class TTracker;

  class CalTimePeakFinder : public art::EDFilter {
  public:

    struct TimePeakHist_t {
      TH1F*  nhits;           // number of hits on a helix  
      TH1F*  energy[2];   
      TH1F*  time[2];
      TH2F*  nhitsvstime;
      TH2F*  nhitsvsenergy;
    };

    struct Hist_t {
      TimePeakHist_t  timePeak;  // helix fit histograms
  
      TH1F* nseeds[2];
    };


  protected:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    //    TStopwatch*   fStopwatch;

    unsigned         _iev;
					// configuration parameters
    int              _diagLevel; 
    int              _debugLevel;
    int              _printfreq;
    int              _useAsFilter; //allows to use the module as a produer or as a filter
//-----------------------------------------------------------------------------
// event object labels
//-----------------------------------------------------------------------------
    std::string      _shLabel ; // MakeStrawHit label (makeSH)
    std::string      _shfLabel;
    std::string      _ccmLabel; // caloClusterModuleLabel

    //    std::string      _dtspecpar;

    StrawHitFlag     _hsel;
    StrawHitFlag     _bkgsel;

    double           _mindt;
    double           _maxdt;
					// time spectrum parameters
    int              _minnhits;

    double           _minClusterEnergy;	// min seed energy
    int              _minClusterSize;   // min size of the seeding cluster
    double           _minClusterTime;   // min time of the seeding cluster

    double           _pitchAngle;
					// outlier cuts
//-----------------------------------------------------------------------------
// cache of event objects
//-----------------------------------------------------------------------------
    art::Handle<mu2e::StrawHitCollection> _strawhitsH;

    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const CaloClusterCollection*          _ccCollection;
    art::Handle<CaloClusterCollection>    _ccH;

    double                                _dtoffset;

    CalTimePeakCollection*                _tpeaks;   // cache of time peaks

    const TTracker*                       _tracker;     // straw tracker geometry
    const Calorimeter*                    _calorimeter; // cached pointer to the calorimeter geometry

    const TrackerCalibrations*            _trackerCalib;

//-----------------------------------------------------------------------------
// diagnostics histograms
//-----------------------------------------------------------------------------
    Hist_t                                _hist;

//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:

    explicit CalTimePeakFinder(const fhicl::ParameterSet& PSet);
    virtual ~CalTimePeakFinder();
    
    virtual void beginJob ();
    virtual bool beginRun (art::Run&);
    virtual bool filter   (art::Event& e);
    virtual void endJob   ();
//-----------------------------------------------------------------------------
// helper functions
//-----------------------------------------------------------------------------
    bool findData         (const art::Event& e);
    void findTimePeaks    (CalTimePeakCollection* TimePeakColl,
			   TimeClusterCollection& OutSeeds);

    void bookHistograms   ();
    void initTimeCluster    (TimeCluster   &TrackSeed   , 
			     CalTimePeak &TPeak       ,
			     int         &ClusterIndex);
    
  };
}
#endif

