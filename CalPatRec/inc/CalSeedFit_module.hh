///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CalPatRec_CalSeedFit_module
#define CalPatRec_CalSeedFit_module

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
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/HelixSeedCollection.hh"
#include "RecoDataProducts/inc/KalSeed.hh"

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
//#include "CalPatRec/inc/CalTimePeak.hh"
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
#include "DataProducts/inc/Helicity.hh"


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

  class CalSeedFit : public art::EDFilter {
  public:
    struct SeedFitHist_t {
      TH1F*  seeddoca[3];
      TH1F*  nhits;           // number of hits on a htrack candidate
      TH1F*  chi2[2];
      TH1F*  p [2];
      TH1F*  NpointsSeed   [2]; //

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
    int              _useAsFilter; //allows to use the module as a produer or as a filter
    int              _rescueHits;
//-----------------------------------------------------------------------------
// event object labels
//-----------------------------------------------------------------------------
    std::string      _shLabel ; // MakeStrawHit label (makeSH)
    std::string      _shDigiLabel;
    std::string      _shpLabel;
    std::string      _shfLabel;
    std::string      _helixSeedLabel;

    double           _maxdtmiss;
					// outlier cuts
    double           _maxadddoca;
    TrkParticle      _tpart;	        // particle type being searched for
    TrkFitDirection  _fdir;		// fit direction in search
    std::vector<double>                   _perr; // diagonal parameter errors to use in the fit


    int              _nhits_from_gen;
//-----------------------------------------------------------------------------
// cache of event objects
//-----------------------------------------------------------------------------
    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;
    const PtrStepPointMCVectorCollection* _listOfMCStrawHits;

    const HelixSeedCollection*            _helixSeeds;

    art::Handle<HelixSeedCollection>      _helixSeedsHandle;

    KalFitHack                            _seedfit;  // Kalman filter config for the Seed fit ( fit using hit wires)

    KalFitResult*                         _sfresult; // seed fit result

    std::vector<StrawHitIndex>            _hitIndices, _goodhits;
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
    Helicity                              _helicity; // cached value of helicity expected for this fit
    double                                _amsign;   // cached sign of angular momentum WRT the z axis 
    CLHEP::HepSymMatrix                   _hcovar; // cache of parameter error covariance matrix

//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    enum fitType {helixFit=0,seedFit,kalFit};
    explicit CalSeedFit(const fhicl::ParameterSet& PSet);
    virtual ~CalSeedFit();
    
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

    void fillSeedFitHistograms(KalFitResult& SFResult);

    void init             (KalFitResult*&  KRes, TrkDefHack* TDef);

  };
}
#endif

