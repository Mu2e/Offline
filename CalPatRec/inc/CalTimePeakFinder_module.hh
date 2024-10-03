///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CalPatRec_CalTimePeakFinder_module
#define CalPatRec_CalTimePeakFinder_module

#ifdef __GCCXML__A
namespace art {
  class EDProducer;
  class Run;
  class Event;
};
#else
#  include "art/Framework/Core/EDProducer.h"
#  include "art/Framework/Principal/Event.h"
#endif

// data
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/HelixVal.hh"

#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawHit.hh"

#include "Offline/RecoDataProducts/inc/KalRepCollection.hh"
#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "Offline/BTrkData/inc/Doublet.hh"

// BaBar
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BField/BField.hh"
#include "BTrk/BField/BFieldFixed.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "Offline/BFieldGeom/inc/BFieldConfig.hh"
#include "BTrk/BaBar/BbrStringUtils.hh"
#include "Offline/BTrkData/inc/TrkStrawHit.hh"
#include "BTrk/TrkBase/HelixParams.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "Offline/TrkPatRec/inc/TrkHitFilter.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"

// #include "CalPatRec/inc/CalTimePeak.hh"
#include "Offline/CalPatRec/inc/CalTimePeakFinder_types.hh"

// Mu2e

//CLHEP
#include "CLHEP/Units/PhysicalConstants.h"

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
  using namespace CalTimePeakFinderTypes;

  class Calorimeter;
  class Tracker;
  class ModuleHistToolBase;

  class CalTimePeakFinder: public art::EDProducer {
  protected:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    unsigned         _iev;
                                        // configuration parameters
    int              _diagLevel;
    int              _debugLevel;
    int              _printfreq;
//-----------------------------------------------------------------------------
// event object labels
//-----------------------------------------------------------------------------
    std::string      _shLabel ; // MakeStrawHit label (makeSH)
    std::string      _shfLabel;
    // std::string      _shpLabel;
    std::string      _ccmLabel; // caloClusterModuleLabel

    StrawHitFlag     _hsel;
    StrawHitFlag     _bkgsel;

    double           _mindt;
    double           _maxdt;
                                        // time spectrum parameters
    int              _minNHits;

    double           _minClusterEnergy;        // min seed energy
    int              _minClusterSize;   // min size of the seeding cluster

    double           _pitchAngle;
    double           _sinPitch;
    double           _beta;
                                        // outlier cuts
//-----------------------------------------------------------------------------
// cache of event objects
//-----------------------------------------------------------------------------
    art::Handle<CaloClusterCollection>    _ccH; // data member, as used from different places

    // const StrawHitCollection*             _shcol;
    // const StrawHitFlagCollection*         _shfcol;
    // const CaloClusterCollection*          _ccCollection;

    double                                _dtoffset;

    const Tracker*                        _tracker;     // straw tracker geometry
    const Calorimeter*                    _calorimeter; // cached pointer to the calorimeter geometry

    const CaloCluster*                     cl;
//-----------------------------------------------------------------------------
// diagnostics
//-----------------------------------------------------------------------------
    Data_t                                 _data;
    std::unique_ptr<ModuleHistToolBase>    _hmanager;
//-----------------------------------------------------------------------------
// functions
//-----------------------------------------------------------------------------
  public:

    explicit CalTimePeakFinder(const fhicl::ParameterSet& PSet);
    virtual ~CalTimePeakFinder();

    virtual void beginJob ();
    virtual void beginRun (art::Run&);
    virtual void produce  (art::Event& e);
    virtual void endJob   ();
//-----------------------------------------------------------------------------
// helper functions
//-----------------------------------------------------------------------------
    bool findData         (const art::Event& e);
    void findTimePeaks    (//CalTimePeakCollection* TimePeakColl,
                           TimeClusterCollection& OutSeeds);

    // void initTimeCluster  (TimeCluster &TrackSeed   ,
    //                            CalTimePeak &TPeak       ,
    //                            int         &ClusterIndex);
  };
}
#endif

