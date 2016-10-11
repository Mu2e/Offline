///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CalPatRec_CalTrkFit_module
#define CalPatRec_CalTrkFit_module

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

#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/TrackSeed.hh"
#include "RecoDataProducts/inc/TrackSeedCollection.hh"


#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruth.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"

// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BaBar/BbrStringUtils.hh"
#include "CalPatRec/inc/TrkDefHack.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkPoca.hh"


#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "TrkPatRec/inc/TrkHitFilter.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"
#include "CalPatRec/inc/CalTimePeak.hh"
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

  class CalTrkFit : public art::EDFilter {
  public:
    struct SeedFitHist_t {
      TH1F*  nhits;           // number of hits on a htrack candidate
      TH1F*  chi2[2];
      TH1F*  p [2];
    };

    struct Hist_t {
      SeedFitHist_t  seedFit;  // helix fit histograms

      TH1F*          ntracks[2];
    };

    Ref*    _ref;

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
    bool             _addhits; 
    bool             _doKalmanFit;
//-----------------------------------------------------------------------------
// event object labels
//-----------------------------------------------------------------------------
    std::string      _shLabel ; // MakeStrawHit label (makeSH)
    std::string      _shpLabel;
    std::string      _shfLabel;
    std::string      _trkseedLabel;
    std::string      _tpeaksLabel;

    double           _maxdtmiss;
					// outlier cuts
    double           _maxadddoca;
    double           _maxaddchi;
    TrkParticle      _tpart;	        // particle type being searched for
    TrkFitDirection  _fdir;		// fit direction in search

    int              _nhits_from_gen;
//-----------------------------------------------------------------------------
// cache of event objects
//-----------------------------------------------------------------------------
    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;
    const PtrStepPointMCVectorCollection* _listOfMCStrawHits;

    const TrackSeedCollection*            _trkseeds;
    const CalTimePeakCollection*          _tpeaks;

    KalFitHack                            _seedfit;  // Kalman filter config for the Seed fit ( fit using hit wires)
    KalFitHack                            _kfit;     // full-blown src/Kalman filter

    KalFitResult*                         _sfresult; // seed fit result
    KalFitResult*                         _kfresult; // full fit result

    std::string                           _iname;	// data instance name

    std::string                           _iname_seed;// data instance name for the output 
					              // of seed fit (used for diagnostics only)

    std::vector<StrawHitIndex>                 _hitIndices;
    int                                   _nindex;
    int                                   _nrescued;    // by the seed fit

    const TTracker*                       _tracker;     // straw tracker geometry
    const Calorimeter*                    _calorimeter; // cached pointer to the calorimeter geometry

    const TrackerCalibrations*            _trackerCalib;

    TFolder*                              _folder;
    int                                   _eventid;
    int                                   _ntracks[2];
//-----------------------------------------------------------------------------
// diagnostics histograms
//-----------------------------------------------------------------------------
    Hist_t                                _hist;

    THackData*                            fHackData;

    int                                   _minNMCHits;

    double                                _mbtime;               // period of 1 microbunch
    SimParticleTimeOffset*                fgTimeOffsets;

    HelixTraj*                            _helTraj;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    enum fitType {helixFit=0,seedFit,kalFit};
    explicit CalTrkFit(const fhicl::ParameterSet& PSet);
    virtual ~CalTrkFit();
    
    virtual void beginJob();
    virtual bool beginRun(art::Run&);
    virtual bool filter (art::Event& event ); 
    virtual void endJob();
//-----------------------------------------------------------------------------
// helper functions
//-----------------------------------------------------------------------------
    bool findData         (const art::Event& e);
//----------------------------------------------------------------------
// 2015 - 02 - 16 Gianipez added the two following functions
//----------------------------------------------------------------------
    void findDoublets     (KalRep* krep, DoubletCollection *dcol);//search doublets in a giventimepeak
    void findLoopApex     (){}//search the straw hits src/closer to the apexes of the helix loops

    void findMissingHits  (KalFitResult& kalfit, std::vector<StrawHitIndex>& indices);
    void bookHistograms   ();

    void init             (KalFitResult*&  KRes, TrkDefHack* TDef);

  };
}
#endif

