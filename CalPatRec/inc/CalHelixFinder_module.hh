///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CalPatRec_CalHelixFinder_module
#define CalPatRec_CalHelixFinder_module

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Services/Optional/TFileService.h"

// data
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/HelixVal.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"

#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BField/BField.hh"
#include "BTrk/BField/BFieldFixed.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BaBar/BbrStringUtils.hh"
#include "BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"

#include "Mu2eBTrk/inc/BaBarMu2eField.hh"
#include "BFieldGeom/inc/BFieldConfig.hh"
#include "TrkPatRec/inc/StrawHitInfo.hh"

#include "CalPatRec/inc/CalHelixFinder_types.hh"
#include "CalPatRec/inc/CalHelixFinderAlg.hh"

// Mu2e
#include "ConditionsService/inc/TrackerCalibrations.hh"

//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"
// root 
#include "TROOT.h"
#include "TFolder.h"
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TApplication.h"
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

  class Calorimeter;
  class TTracker;
  class CprModuleHistBase;

  class CalHelixFinder : public art::EDFilter {
  protected:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
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
    art::Handle<StrawHitCollection>       _strawhitsH;
    art::Handle<TimeClusterCollection>    _timeclcolH;

    fhicl::ParameterSet*                  _timeOffsets;

    const StrawHitCollection*             _shcol;
    const StrawHitFlagCollection*         _shfcol;
    const StrawHitPositionCollection*     _shpcol;
    const TimeClusterCollection*          _timeclcol;

    HelixTraj*                            _helTraj;
    CalHelixFinderAlg                     _hfinder;	

    const TTracker*                       _tracker     ; // straw tracker geometry
    const TrackerCalibrations*            _trackerCalib; // 
    const Calorimeter*                    _calorimeter ; // cached pointer to the calorimeter geometry
//-----------------------------------------------------------------------------
// diagnostics 
//-----------------------------------------------------------------------------
    CalHelixFinder_Hist_t                 _hist;
    CalHelixFinder_Data_t                 _data;

    std::unique_ptr<CprModuleHistBase>    _hmanager;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:
    enum fitType {helixFit=0,seedFit,kalFit};
    explicit CalHelixFinder(const fhicl::ParameterSet& PSet);
    virtual ~CalHelixFinder();
    
    virtual void beginJob();
    virtual bool beginRun(art::Run&   run   );
    virtual bool filter  (art::Event& event ); 
    virtual void endJob();
//-----------------------------------------------------------------------------
// helper functions
//-----------------------------------------------------------------------------
    bool findData         (const art::Event& e);
//----------------------------------------------------------------------
// 2015 - 02 - 16 Gianipez added the two following functions
//----------------------------------------------------------------------
    void initHelixSeed      (HelixSeed &TrackSeed, CalHelixFinderData &HfResult);

    int  initHelixFinderData(CalHelixFinderData&                Data,
			     const TrkParticle&                 TPart,
			     const TrkFitDirection&             FDir,
			     const StrawHitCollection*          StrawCollection ,
			     const StrawHitPositionCollection*  ShPosCollection , 
			     const StrawHitFlagCollection*      ShFlagCollection);

  };
}
#endif

