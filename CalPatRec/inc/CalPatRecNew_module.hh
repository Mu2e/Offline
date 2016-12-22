///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CalPatRec_CalPatRecNew_module
#define CalPatRec_CalPatRecNew_module

#ifdef __GCCXML__A
namespace art {
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
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/HelixSeedCollection.hh"



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

  class CalPatRecNew : public art::EDFilter {
  public:
    struct HelixFitHist_t {
      TH1F*  nhits;           // number of hits on a helix  
      TH1F*  radius[2];   
      TH1F*  chi2XY[2];
      TH1F*  chi2ZPhi[2];
      TH1F*  pT[2];
      TH1F*  p [2];
      TH2F*  nhitsvspT;
      TH2F*  nhitsvsp;
    };

    struct Hist_t {
      HelixFitHist_t  helixFit;  // helix fit histograms

      TH1F* nseeds[2];
    };

    Ref*                                   _ref;

  protected:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    //    TStopwatch*   fStopwatch;

    unsigned                              _iev;
					// configuration parameters
    int                                   _diagLevel; 
    int                                   _debugLevel;
    int                                   _printfreq;
    int                                   _useAsFilter; //allows to use the module as a produer or as a filter
//-----------------------------------------------------------------------------
// event object labels
//-----------------------------------------------------------------------------
    std::string                           _shLabel ; // MakeStrawHit label (makeSH)
    std::string                           _shpLabel;
    std::string                           _shfLabel;
    std::string                           _timeclLabel;

    TrkParticle                           _tpart;	        // particle type being searched for
    TrkFitDirection                       _fdir;		// fit direction in search

//-----------------------------------------------------------------------------
// cache of event objects
//-----------------------------------------------------------------------------
    art::Handle<mu2e::StrawHitCollection>    _strawhitsH;
    art::Handle<mu2e::TimeClusterCollection> _timeclcolH;

    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;
    const TimeClusterCollection*          _timeclcol;
    const CalTimePeakCollection*          _tpeaks;

    HelixFitHack                          _hfit;	

    XYZPHackVector                        _index;
    int                                   _nindex;

    const TTracker*                       _tracker;     // straw tracker geometry
    const Calorimeter*                    _calorimeter; // cached pointer to the calorimeter geometry

    const TrackerCalibrations*            _trackerCalib;

    TFolder*                              _folder;
//-----------------------------------------------------------------------------
// diagnostics histograms
//-----------------------------------------------------------------------------
    Hist_t                                _hist;

    THackData*                            fHackData;

    //    double                   _mbtime;               // period of 1 microbunch
    //    SimParticleTimeOffset*   fgTimeOffsets;

    HelixTraj*                            _helTraj;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    enum fitType {helixFit=0,seedFit,kalFit};
    explicit CalPatRecNew(const fhicl::ParameterSet& PSet);
    virtual ~CalPatRecNew();
    
    virtual void beginJob();
    virtual bool beginRun(art::Run&);
    virtual bool filter  (art::Event& event ); 
    virtual void endJob();
//-----------------------------------------------------------------------------
// helper functions
//-----------------------------------------------------------------------------
    bool findData         (const art::Event& e);

//----------------------------------------------------------------------
// 2015 - 02 - 16 Gianipez added the two following functions
//----------------------------------------------------------------------

    void bookHistograms   ();
    void initHelixSeed    (HelixSeed                             &TrackSeed , 
			   TrkDefHack                            &Seeddef   ,
			   HelixFitHackResult                    &HfResult  ,
			   art::Ptr<TimeCluster>                 &TPeak     );

  };
}
#endif

