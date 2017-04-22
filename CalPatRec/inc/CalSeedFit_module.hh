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
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
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
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkPoca.hh"


#include "RecoDataProducts/inc/KalRepCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"

#include "CalPatRec/inc/KalFitHackNew.hh"
#include "CalPatRec/inc/KalFitResultNew.hh"
#include "CalPatRec/inc/CalSeedFit_types.hh"
#include "CalPatRec/inc/CprModuleHistBase.hh"

// Mu2e
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

namespace fhicl {
  class ParameterSet;
}

namespace mu2e {

  class Calorimeter;
  class TTracker;

  class CalSeedFit : public art::EDFilter {
  protected:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    unsigned            _iev;
					// configuration parameters
    int                 _diagLevel; 
    int                 _debugLevel;
    int                 _printfreq;
    int                 _useAsFilter;   // 0: producer, 1: filter
    int                 _rescueHits;
//-----------------------------------------------------------------------------
// event object labels
//-----------------------------------------------------------------------------
    std::string         _shLabel ;      // MakeStrawHit label (makeSH)
    std::string         _shDigiLabel;
    std::string         _shpLabel;
    std::string         _shfLabel;
    std::string         _helixSeedLabel;

    double              _maxdtmiss;
					// outlier cuts
    double              _maxAddDoca;
    double              _maxAddChi;
    TrkParticle         _tpart;	        // particle type being searched for
    TrkFitDirection     _fdir;		// fit direction in search
    std::vector<double> _perr;          // diagonal parameter errors to use in the fit

    int                 _nhits_from_gen;//
//-----------------------------------------------------------------------------
// cache of event objects
//-----------------------------------------------------------------------------
    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;
    int                                   _nhits;  // N hits in _shcol

    const PtrStepPointMCVectorCollection* _listOfMCStrawHits;

    const HelixSeedCollection*            _helixSeeds;

    art::Handle<HelixSeedCollection>      _helixSeedsHandle;

    KalFitHackNew                         _fitter;  // Kalman filter config for the Seed fit ( fit using hit wires)

    KalFitResultNew                       _result; // seed fit result

    const TTracker*                       _tracker;     // straw tracker
    const TrackerCalibrations*            _trackerCalib;
    const Calorimeter*                    _calorimeter; // cached pointer to the calorimeter

    TFolder*                              _folder;
    int                                   _eventid;
    int                                   _ntracks[2];
//-----------------------------------------------------------------------------
// diagnostics 
//-----------------------------------------------------------------------------
    CalSeedFit_Hist_t                     _hist;
    CalSeedFit_Data_t                     _data;
    std::unique_ptr<CprModuleHistBase>    _hmanager;

    double                                _amsign;   // cached sign of angular momentum WRT the z axis 
    CLHEP::HepSymMatrix                   _hcovar;   // cache of parameter error covariance matrix
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    enum fitType { helixFit=0, seedFit=1, kalFit=2 };
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
    void findMissingHits(KalFitResultNew& KRes);
  };
}
#endif

