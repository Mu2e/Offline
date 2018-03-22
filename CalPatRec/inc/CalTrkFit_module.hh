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
#include "art/Framework/Services/Optional/TFileService.h"
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
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHit.hh"
#include "RecoDataProducts/inc/StrawHitIndex.hh"
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
//#include "CalPatRec/inc/CalTimePeak.hh"
#include "RecoDataProducts/inc/Doublet.hh"

#include "TrkReco/inc/DoubletAmbigResolver.hh"

#include "TROOT.h"
#include "TFolder.h"
#include "CalPatRec/inc/KalFitHackNew.hh"
#include "CalPatRec/inc/CalTrkFit_types.hh"
#include "CalPatRec/inc/ModuleHistToolBase.hh"

//#include "CalPatRec/inc/THackData.hh"

// Mu2e
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

namespace fhicl {
  class ParameterSet;
}

namespace mu2e {
  using namespace CalTrkFitTypes;
  
  class Calorimeter;
  class TTracker;

  class CalTrkFit : public art::EDFilter {
  protected:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    unsigned            _iev;
	    	  	  		// configuration parameters
    int                 _diagLevel; 
    int                 _debugLevel;
    int                 _printfreq;
    int                 _useAsFilter;   // allows to use the module as a produer or as a filter
    bool                _addhits; 
    std::vector<double> _zsave;
//-----------------------------------------------------------------------------
// event object labels
//-----------------------------------------------------------------------------
    std::string         _shLabel ;      // MakeStrawHit label (makeSH)
    std::string         _shDigiLabel;
    std::string         _shpLabel;
    std::string         _shfLabel;
    std::string         _trkseedLabel;
    std::string         _tpeaksLabel;

    double              _maxdtmiss;
					// outlier cuts
    double              _maxadddoca;
    double              _maxaddchi;
    TrkFitFlag          _goodseed;
    
    TrkParticle         _tpart;	        // particle type being searched for
    TrkFitDirection     _fdir;		// fit direction in search
//-----------------------------------------------------------------------------
// cache of event objects
//-----------------------------------------------------------------------------
    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;

    const KalSeedCollection*              _trkseeds;

    KalFitHackNew                         _fitter;      // full-blown src/Kalman filter

    KalFitResultNew                       _result;      // full fit result

    std::vector<StrawHitIndex>            _hitIndices;
    int                                   _nindex;
    int                                   _nrescued;    // by the seed fit

    const TrackerCalibrations*            _trackerCalib;

    TFolder*                              _folder;
    int                                   _eventid;
//-----------------------------------------------------------------------------
// diagnostics histograms
//-----------------------------------------------------------------------------
    Data_t                                _data;
    std::unique_ptr<ModuleHistToolBase>   _hmanager;
    vector<Doublet>*                      _listOfDoublets;

    double                                _mbtime;      // period of 1 microbunch
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
    void findMissingHits  (KalFitResultNew&  KRes);
  };
}
#endif

