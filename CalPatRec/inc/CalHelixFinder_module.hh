///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CalPatRec_CalHelixFinder_module
#define CalPatRec_CalHelixFinder_module

#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"

// data
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloHit.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/HelixVal.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"

#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/StereoHit.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawHit.hh"

#include "DataProducts/inc/Helicity.hh"

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

#include "BFieldGeom/inc/BFieldConfig.hh"

#include "CalPatRec/inc/CalHelixFinder_types.hh"
#include "CalPatRec/inc/CalHelixFinderAlg.hh"
#include "CalPatRec/inc/CalHelixFinderData.hh"

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
  using namespace CalHelixFinderTypes;

  class Calorimeter;
  class Tracker;
  class ModuleHistToolBase;

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
    // std::string                           _shpLabel;
    std::string                           _shfLabel;
    std::string                           _timeclLabel;
    
    int                                   _minNHitsTimeCluster; //min nhits within a TimeCluster after check of Delta-ray hits

    TrkParticle                           _tpart;	        // particle type being searched for
    TrkFitDirection                       _fdir;		// fit direction in search
//-----------------------------------------------------------------------------
// cache of event objects
//-----------------------------------------------------------------------------
    art::Handle<ComboHitCollection>       _strawhitsH;
    art::Handle<TimeClusterCollection>    _timeclcolH;

    fhicl::ParameterSet*                  _timeOffsets;

    const ComboHitCollection*             _chcol;
    const StrawHitFlagCollection*         _shfcol;
    // const StrawHitPositionCollection*     _shpcol;
    const TimeClusterCollection*          _timeclcol;

    HelixTraj*                            _helTraj;
    CalHelixFinderAlg                     _hfinder;	
    CalHelixFinderData                    _hfResult;
    std::vector<mu2e::Helicity>           _hels; // helicity values to fit

    double                                _bz0;
    const Tracker*                        _tracker     ; // straw tracker geometry
    const Calorimeter*                    _calorimeter ; // cached pointer to the calorimeter geometry
//-----------------------------------------------------------------------------
// diagnostics 
//-----------------------------------------------------------------------------
    CalHelixFinderTypes::Data_t           _data;

    std::unique_ptr<ModuleHistToolBase>   _hmanager;
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
			     const ComboHitCollection*          ComboCollection ,
			     // const StrawHitPositionCollection*  ShPosCollection , 
			     const StrawHitFlagCollection*      ShFlagCollection);
    
    int  goodHitsTimeCluster(const TimeCluster* TimeCluster);
    
    void pickBestHelix(std::vector<HelixSeed>& HelVec, int &Index_best);
  };
}
#endif

