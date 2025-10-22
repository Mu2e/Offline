///////////////////////////////////////////////////////////////////////////////
//
///////////////////////////////////////////////////////////////////////////////
#ifndef CalPatRec_CalHelixFinder_module
#define CalPatRec_CalHelixFinder_module

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"

// data
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/HelixVal.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"

#include "Offline/RecoDataProducts/inc/StrawHit.hh"

#include "Offline/DataProducts/inc/Helicity.hh"

// BaBar
#include "BTrkLegacy/inc/HelixTraj.hh"
#include "Offline/BFieldGeom/inc/BFieldConfig.hh"

#include "Offline/CalPatRec/inc/CalHelixFinder_types.hh"
#include "Offline/CalPatRec/inc/CalHelixFinderAlg.hh"
#include "Offline/CalPatRec/inc/CalHelixFinderData.hh"

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

  class CalHelixFinder : public art::EDProducer {
  protected:
//-----------------------------------------------------------------------------
// data members
//-----------------------------------------------------------------------------
    unsigned                              _iev;
                                                        // configuration parameters
    int                                   _diagLevel;
    int                                   _debugLevel;
    int                                   _printfreq;
//-----------------------------------------------------------------------------
// event object labels
//-----------------------------------------------------------------------------
    std::string                           _shLabel ; // MakeStrawHit label (makeSH)
    // std::string                           _shpLabel;
    std::string                           _timeclLabel;

    int                                   _minNHitsTimeCluster; //min nhits within a TimeCluster after check of Delta-ray hits

    int                                   _fitparticle;
    TrkParticle                           _tpart;                // particle type being searched for
    TrkFitDirection                       _fdir;                // fit direction in search
    bool                                  _doSingleOutput;
    float                                 _maxEDepAvg;
//-----------------------------------------------------------------------------
// cache of event objects
//-----------------------------------------------------------------------------
    art::Handle<ComboHitCollection>       _strawhitsH;
    art::Handle<TimeClusterCollection>    _timeclcolH;

    fhicl::ParameterSet*                  _timeOffsets;

    const ComboHitCollection*             _chcol;
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


    struct Config
    {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<int>                           diagLevel{            Name("diagLevel"),                  Comment("Diag"),0 };
      fhicl::Atom<int>                           debugLevel{           Name("debugLevel"),                 Comment("Debug"),0 };
      fhicl::Atom<int>                           printfreq{            Name("printFrequency"),                  Comment("Print Frequency") };
      fhicl::Atom<std::string>                   shLabel{              Name("StrawHitCollectionLabel"),                    Comment("StrawHit Collection Label") };
      fhicl::Atom<std::string>                   timeclLabel{          Name("TimeClusterCollectionLabel"),                Comment("TimeCluster Collection Label") };
      fhicl::Atom<int>                           minNHitsTimeCluster{  Name("minNHitsTimeCluster"),        Comment("Min NHits in TimeCluster") };
      fhicl::Atom<int>                           fitparticle{          Name("fitparticle"),                      Comment("Particle Type Searched For") };
      fhicl::Atom<std::string>                   fitdirection{         Name("fitdirection"),               Comment("Fit Direction in Search (\"downstream\" or \"upstream\")") };
      fhicl::Atom<bool>                          doSingleOutput{       Name("doSingleOutput"),             Comment("Do Single Output") };
      fhicl::Atom<float>                         maxEDepAvg{           Name("maxEDepAvg"),                 Comment("Max Avg EDep") };
      fhicl::Table<CalHelixFinderAlg::Config>    hfinder{              Name("HelixFinderAlg"),                    Comment("CalHelixFinderAlg Config") };
      fhicl::Table<CalHelixFinderTypes::Config>  diagPlugin{           Name("diagPlugin"),                 Comment("Diag Plugin") };
      fhicl::Sequence<int>                       Helicities{           Name("Helicities"),                 Comment("Helicity values") };
    };

    enum fitType {helixFit=0,seedFit,kalFit};

    explicit CalHelixFinder(const art::EDProducer::Table<Config>& config);
    virtual ~CalHelixFinder();

    virtual void beginJob();
    virtual void beginRun(art::Run&   run   );
    virtual void produce (art::Event& event );
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
                             const ComboHitCollection*          ComboCollection);

    int  goodHitsTimeCluster(const TimeCluster* TimeCluster);

    void pickBestHelix(std::vector<HelixSeed>& HelVec, int &Index_best);
  };
}
#endif
