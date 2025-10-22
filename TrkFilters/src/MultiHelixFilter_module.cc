//
//  Filter to select multi-helix events
//  Michael MacKenzie, 2025
//
// framework
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "fhiclcpp/types/OptionalSequence.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"

#include "Offline/BFieldGeom/inc/BFieldManager.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

// data
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitIndex.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/DataProducts/inc/Helicity.hh"
// mu2e
#include "Offline/Mu2eUtilities/inc/HelixTool.hh"
// helper function
#include "Offline/GeneralUtilities/inc/PhiPrescalingParams.hh"
#include "Offline/GeneralUtilities/inc/ParameterSetHelpers.hh"

// ROOT
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"

// c++
#include <string>
#include <vector>
#include <set>
#include <limits>
#include <iostream>
#include <memory>

using namespace CLHEP;

namespace mu2e
{
  class MultiHelixFilter : public art::EDFilter
  {
  public:

    //--------------------------------------------------------------------------------------
    struct HelixCutsConfig {
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::OptionalAtom<double>  minMomentum        {Name("minMomentum"),        Comment("minMomentum")                                       };
      fhicl::OptionalAtom<double>  maxMomentum        {Name("maxMomentum"),        Comment("maxMomentum")                                       };
      fhicl::OptionalAtom<double>  minPt              {Name("minPt"),              Comment("minPt")                                             };
      fhicl::OptionalAtom<bool>    requireCaloCluster {Name("requireCaloCluster"), Comment("requireCaloCluster")                                };
      fhicl::OptionalAtom<int>     minNStrawHits      {Name("minNStrawHits"),      Comment("minNStrawHits")                                     };
      fhicl::OptionalAtom<double>  minHitRatio        {Name("minHitRatio"),        Comment("minHitRatio")                                       };
      fhicl::OptionalAtom<double>  maxChi2XY          {Name("maxChi2XY"),          Comment("maxChi2XY")                                         };
      fhicl::OptionalAtom<double>  maxChi2PhiZ        {Name("maxChi2PhiZ"),        Comment("maxChi2PhiZ")                                       };
      fhicl::OptionalAtom<double>  minD0              {Name("minD0"),              Comment("minD0")                                             };
      fhicl::OptionalAtom<double>  maxD0              {Name("maxD0"),              Comment("maxD0")                                             };
      fhicl::OptionalAtom<double>  minAbsLambda       {Name("minAbsLambda"),       Comment("minAbsLambda")                                      };
      fhicl::OptionalAtom<double>  maxAbsLambda       {Name("maxAbsLambda"),       Comment("maxAbsLambda")                                      };
      fhicl::OptionalAtom<double>  minNLoops          {Name("minNLoops"),          Comment("minNLoops")                                         };
      fhicl::OptionalAtom<double>  maxNLoops          {Name("maxNLoops"),          Comment("maxNLoops")                                         };
      fhicl::OptionalAtom<double>  minSlopeSig        {Name("minSlopeSig"),        Comment("Minimum helix seed slope significance")             };
      fhicl::OptionalAtom<double>  maxSlopeSig        {Name("maxSlopeSig"),        Comment("Maximum helix seed slope significance")             };
      fhicl::OptionalAtom<int>     helicity           {Name("helicity"),           Comment("helicity")                                          };
      fhicl::Sequence<std::string> helixFitFlag       {Name("helixFitFlag"),       Comment("helixFitFlag"), std::vector<std::string>{"HelixOK"} };
    };

    //--------------------------------------------------------------------------------------
    struct HelixCuts {

      HelixCuts() {
      }

      HelixCuts(const HelixCutsConfig& config) : _goodh(config.helixFitFlag()) {
        if(config.minMomentum       ()) _minMomentum              = config.minMomentum       ().value();
        if(config.maxMomentum       ()) _maxMomentum              = config.maxMomentum       ().value();
        if(config.minPt             ()) _minPt                    = config.minPt             ().value();
        if(config.requireCaloCluster()) _requireCaloCluster       = config.requireCaloCluster().value();
        if(config.minNStrawHits     ()) _minNStrawHits            = config.minNStrawHits     ().value();
        if(config.minHitRatio       ()) _minHitRatio              = config.minHitRatio       ().value();
        if(config.maxChi2XY         ()) _maxChi2XY                = config.maxChi2XY         ().value();
        if(config.maxChi2PhiZ       ()) _maxChi2PhiZ              = config.maxChi2PhiZ       ().value();
        if(config.minD0             ()) _minD0                    = config.minD0             ().value();
        if(config.maxD0             ()) _maxD0                    = config.maxD0             ().value();
        if(config.minAbsLambda      ()) _minAbsLambda             = config.minAbsLambda      ().value();
        if(config.maxAbsLambda      ()) _maxAbsLambda             = config.maxAbsLambda      ().value();
        if(config.minNLoops         ()) _minNLoops                = config.minNLoops         ().value();
        if(config.maxNLoops         ()) _maxNLoops                = config.maxNLoops         ().value();
        if(config.minSlopeSig       ()) _minSlopeSig              = config.minSlopeSig       ().value();
        if(config.maxSlopeSig       ()) _maxSlopeSig              = config.maxSlopeSig       ().value();
        if(config.helicity          ()) _hel                      = config.helicity          ().value();
      }
      void setB(const float bz0) { _bz0 = bz0; }
      void setTracker(const Tracker* tracker) { _tracker = tracker; }
      bool applyCuts(const HelixSeed* helix) {
        if(_debugLevel > 2) std::cout << "Testing helix:\n";

        // helix exists and has a good status
        if(!helix) return false;
        if(!helix->status().hasAllProperties(_goodh))                return false;

        // momentum selections
        constexpr float mmToMeV = 3./10.;
        const float mom = _bz0 * mmToMeV * helix->helix().momentum(); //convert mm --> momentum using B-field
        const float pt  = _bz0 * mmToMeV * helix->helix().radius  (); //convert mm --> momentum using B-field
        if(_debugLevel > 2) std::cout << "  p = " << mom << " pT = " << pt << std::endl;

        if(_minMomentum > 0. && mom < _minMomentum)                  return false;
        if(_maxMomentum > 0. && mom > _maxMomentum)                  return false;
        if(_minPt       > 0. && pt  < _minPt)                        return false;

        // initialize the helix tool
        HelixTool hTool(helix, _tracker);

        // geometric selections
        const Helicity hel         = helix->helix().helicity();
        const float d0             = helix->helix().rcent() - helix->helix().radius();
        const float lambda         = std::fabs(helix->helix().lambda());
        const float nLoops         = hTool.nLoops();
        const float slope          = helix->recoDir().slope();
        const float slopeErr       = std::fabs(helix->recoDir().slopeErr());
        const float slopeSignedSig = (slopeErr > 0.f) ? slope/slopeErr : 0.f; //signifance from 0 signed by the slope direction
        if(_debugLevel > 2) std::cout << "  hel = " << hel.value() << " d0 = " << d0 << " lambda = " << lambda << " nLoops = " << nLoops
                                      << " slopeSig = " << slopeSignedSig << std::endl;

        if(_hel != 0 && hel == Helicity(_hel))                       return false;
        if(d0 < _minD0)                                              return false;
        if(d0 > _maxD0)                                              return false;
        if(lambda < _minAbsLambda)                                   return false;
        if(lambda > _maxAbsLambda)                                   return false;
        if(nLoops < _minNLoops)                                      return false;
        if(nLoops > _maxNLoops)                                      return false;
        if(slopeSignedSig < _minSlopeSig)                            return false;
        if(slopeSignedSig > _maxSlopeSig)                            return false;

        // quality selections
        const float chi2XY         = helix->helix().chi2dXY();
        const float chi2PhiZ       = helix->helix().chi2dZPhi();
        const int   nhits          = hTool.nstrawhits();
        const float hRatio         = hTool.hitRatio();
        const bool  hasCluster     = helix->caloCluster().isNonnull();
        if(_debugLevel > 2) std::cout << " chi2XY = " << chi2XY << " chi2PhiZ = " << chi2PhiZ
                                      << " nhits = " << nhits << " hitRatio = " << hRatio << " hasCC = " << hasCluster << std::endl;

        if(chi2XY > _maxChi2XY)                                      return false;
        if(chi2PhiZ > _maxChi2PhiZ)                                  return false;
        if(nhits < _minNStrawHits)                                   return false;
        if(hRatio < _minHitRatio)                                    return false;
        if(_requireCaloCluster && !hasCluster)                       return false;

        if(_debugLevel > 2) std::cout << "  --> passed all cuts\n";
        // passes all cuts
        return true;
      }

      double _minMomentum              = -1.   ;
      double _maxMomentum              = -1.   ;
      double _minPt                    = -1.   ;
      bool   _requireCaloCluster       = false ;
      int    _minNStrawHits            = -1.   ;
      double _minHitRatio              = -1.   ;
      double _maxChi2XY                = -1.   ;
      double _maxChi2PhiZ              = -1.   ;
      double _minD0                    = std::numeric_limits<double>::lowest();
      double _maxD0                    = std::numeric_limits<double>::max();
      double _minAbsLambda             = -1.   ;
      double _maxAbsLambda             = -1.   ;
      double _minNLoops                = -1.   ;
      double _maxNLoops                = -1.   ;
      double _minSlopeSig              = std::numeric_limits<double>::lowest();
      double _maxSlopeSig              = std::numeric_limits<double>::max();
      int    _hel                      =  0    ;

      TrkFitFlag _goodh; // helix fit flag

      double _bz0                      =  1.   ;
      const Tracker* _tracker          = nullptr;
      int    _debugLevel               =  0    ;
    };

    //--------------------------------------------------------------------------------------
    struct Config{
      using Name    = fhicl::Name;
      using Comment = fhicl::Comment;
      fhicl::Atom<art::InputTag>                     helixSeedCollection  { Name("helixSeedCollection"), Comment("helixSeedCollection") };
      fhicl::Atom<int>                               debugLevel           { Name("debugLevel"),          Comment("debugLevel")     , 0 };
      fhicl::Atom<bool>                              noFilter             { Name("noFilter"),            Comment("don't apply the filter decision"), false};
      fhicl::Atom<unsigned>                          minNHelices          { Name("minNHelices"),         Comment("minimum number of helices passing the cuts"), 1};
      fhicl::Atom<float>                             timeWindow           { Name("timeWindow"),          Comment("Maximum time between earliest and latest helix"), -1.};
      fhicl::Atom<float>                             maxOverlap           { Name("maxOverlap"),          Comment("Maximum hit overlap between helices"), -1.};
      fhicl::Atom<float>                             minArcGap            { Name("minArcGap"),           Comment("Minimum distance between helix arcs of same curve direction"), -1.};
      fhicl::Sequence<fhicl::Table<HelixCutsConfig>> HelixCuts            { Name("HelixCuts"),           Comment("Helix cut sets")};
      fhicl::Atom<bool>                              doHistograms         { Name("doHistograms"),        Comment("Make debug histograms"), false};
    };


    using Parameters = art::EDFilter::Table<Config>;

    explicit MultiHelixFilter(const Parameters& conf);
    virtual bool filter(art::Event& event) override;
    virtual void beginJob();
    virtual bool beginRun(art::Run&   run   );
    virtual bool endRun( art::Run& run ) override;

  private:

    struct Hist_t {
      TH1* p            ;
      TH1* pt           ;
      TH1* t            ;
      TH1* d0           ;
      TH1* lambda       ;
      TH1* clusterEnergy;
      TH1* nStrawHits   ;
      TH1* hitRatio     ;
      TH1* chi2XY       ;
      TH1* chi2PhiZ     ;
      TH1* nLoops       ;
      TH1* slopeSig     ;
      TH1* helicity     ;

      // multi-helix info
      TH1* overlap      ;
      TH1* arcGap       ;
      TH1* deltaT       ;

      // Per-fill info
      TH1* helicitySum  ;
      TH1* nHelices     ;
      TH1* nInSet       ;
    };

    float evaluateOverlap(const HelixSeed* h1, const HelixSeed* h2);
    float evaluateArcGap(const HelixSeed* h1, const HelixSeed* h2);
    bool checkHelixList(std::vector<const HelixSeed*> helices);
    void doPermutations(const std::vector<size_t>& given, const HelixSeedCollection* hcol, std::set<size_t>& accepted);
    void bookHistograms(const int index, const char* name, const char* title);
    void fillHistograms(Hist_t* Hist, const std::set<size_t>& indices, const HelixSeedCollection* hcol);

    art::InputTag  _hsTag;
    int            _debugLevel;
    bool           _noFilter;
    unsigned       _minNHelices;
    float          _timeWindow;
    float          _maxOverlap;
    float          _minArcGap;
    // cut sets
    std::vector<HelixCuts> _helixCuts;
    bool          _doHistograms;

    double         _bz0;
    const Tracker* _tracker;

    // counters
    unsigned      _nevt, _npass;

    // histogram info
    Hist_t* _hists[10];
  };

  //--------------------------------------------------------------------------------------
  MultiHelixFilter::MultiHelixFilter(const Parameters& config):
    art::EDFilter{config}
    , _hsTag             (config().helixSeedCollection())
    , _debugLevel        (config().debugLevel())
    , _noFilter          (config().noFilter())
    , _minNHelices       (config().minNHelices())
    , _timeWindow        (config().timeWindow())
    , _maxOverlap        (config().maxOverlap())
    , _minArcGap         (config().minArcGap())
    , _nevt              (0)
    , _npass             (0)
    {
      produces<TriggerInfo>();
      for(auto& cuts : config().HelixCuts()) _helixCuts.push_back(cuts);
    }

  //--------------------------------------------------------------------------------------
  void MultiHelixFilter::beginJob() {
    if(_doHistograms) {
      art::ServiceHandle<art::TFileService> tfs;
      bookHistograms(0, "all_events", "All events");
      bookHistograms(1, "accepted_events", "Accepted events");
    }
  }

  //--------------------------------------------------------------------------------------
  bool MultiHelixFilter::beginRun(art::Run & run){
    // get the B-field z component at the tracker center
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    _bz0 = bfmgr->getBField(vpoint_mu2e).z();

    // get the tracker info
    mu2e::GeomHandle<mu2e::Tracker> th;
    _tracker = th.get();

    // initialize cut parameters
    for(auto& cuts : _helixCuts) {
      cuts.setB(_bz0);
      cuts.setTracker(_tracker);
      cuts._debugLevel = _debugLevel;
    }
    return true;
  }

  //--------------------------------------------------------------------------------------
  void MultiHelixFilter::bookHistograms(const int index, const char* name, const char* title) {
    _hists[index] = new Hist_t;
    auto Hist = _hists[index];
    art::ServiceHandle<art::TFileService> tfs;
    art::TFileDirectory tfdir = tfs->mkdir(name, title);

    Hist->p              = tfdir.make<TH1F>("p"             , "Helix momentum"           , 100,    0.,  200.);
    Hist->pt             = tfdir.make<TH1F>("pt"            , "Helix transverse momentum", 100,    0.,  200.);
    Hist->t              = tfdir.make<TH1F>("t"             , "Helix time"               , 100,    0., 2000.);
    Hist->d0             = tfdir.make<TH1F>("d0"            , "Helix d_{0}"              , 100, -500.,  500.);
    Hist->lambda         = tfdir.make<TH1F>("lambda"        , "Helix |#lambda|"          , 100,    0., 2000.);
    Hist->clusterEnergy  = tfdir.make<TH1F>("clusterEnergy" , "Helix cluster energy"     , 100,    0.,  200.);
    Hist->nStrawHits     = tfdir.make<TH1F>("nStrawHits"    , "Helix N(straw hits)"      , 100,    0.,  100.);
    Hist->hitRatio       = tfdir.make<TH1F>("hitRatio"      , "Helix hit ratio"          , 100,    0.,    2.);
    Hist->chi2XY         = tfdir.make<TH1F>("chi2XY"        , "Helix #chi^{2}(x,y)"      , 100,    0.,   10.);
    Hist->chi2PhiZ       = tfdir.make<TH1F>("chi2PhiZ"      , "Helix #chi^{2}(#phi,z)"   , 100,    0.,   10.);
    Hist->nLoops         = tfdir.make<TH1F>("nLoops"        , "Helix N(loops)"           , 100,    0.,   10.);
    Hist->slopeSig       = tfdir.make<TH1F>("slopeSig"      , "Helix slope/uncertainty"  , 100,  -20.,   20.);
    Hist->helicity       = tfdir.make<TH1D>("helicity"      , "Helix helicity"           ,   4,   -2.,    2.);

    Hist->overlap        = tfdir.make<TH1F>("overlap"       , "Helix overlap"            , 100,    0.,    1.);
    Hist->arcGap         = tfdir.make<TH1F>("arcGap"        , "Helix arc-gap"            , 100,    0.,  500.);
    Hist->deltaT         = tfdir.make<TH1F>("deltaT"        , "Helix #Deltat"            , 100,    0.,  500.);

    Hist->helicitySum    = tfdir.make<TH1F>("helicitySum"   , "#Sigma(helicity)"         ,  20,  -10.,   10.);
    Hist->nHelices       = tfdir.make<TH1F>("nHelices"      , "N(helices)"               ,  20,    0.,   20.);
    Hist->nInSet         = tfdir.make<TH1F>("nInSet"        , "N(helices) in set"        ,  20,    0.,   20.);
  }

  //--------------------------------------------------------------------------------------
  void MultiHelixFilter::fillHistograms(Hist_t* Hist, const std::set<size_t>& indices, const HelixSeedCollection* hcol) {
    std::set<size_t> checked; // to avoid double counting in comparisons
    float helSum = 0.;
    for(auto index : indices) {
      const HelixSeed& helix = hcol->at(index);
      HelixTool hTool(&helix, _tracker);
      constexpr float mmToMeV = 3./10.;

      Hist->p              ->Fill(_bz0 * mmToMeV * helix.helix().momentum());
      Hist->pt             ->Fill(_bz0 * mmToMeV * helix.helix().radius  ());
      Hist->t              ->Fill(helix.t0().t0());
      Hist->d0             ->Fill(helix.helix().rcent() - helix.helix().radius());
      Hist->lambda         ->Fill(std::fabs(helix.helix().lambda()));
      Hist->clusterEnergy  ->Fill((helix.caloCluster().isNonnull()) ? helix.caloCluster()->energyDep() : 0.);
      Hist->nStrawHits     ->Fill(hTool.nstrawhits());
      Hist->hitRatio       ->Fill(hTool.hitRatio());
      Hist->chi2XY         ->Fill(helix.helix().chi2dXY());
      Hist->chi2PhiZ       ->Fill(helix.helix().chi2dZPhi());
      Hist->nLoops         ->Fill(hTool.nLoops());
      Hist->helicity       ->Fill(helix.helix().helicity().value());
      const float slope          = helix.recoDir().slope();
      const float slopeErr       = std::fabs(helix.recoDir().slopeErr());
      const float slopeSignedSig = (slopeErr > 0.f) ? slope/slopeErr : 0.f; //signifance from 0 signed by the slope direction
      Hist->slopeSig       ->Fill(slopeSignedSig);
      helSum += helix.helix().helicity().value();

      for(auto other : indices) {
        if(other == index) continue;
        if(checked.find(other) != checked.end()) continue;
        const HelixSeed& other_helix = hcol->at(other);
        Hist->overlap ->Fill(evaluateOverlap(&helix, &other_helix));
        Hist->arcGap  ->Fill(evaluateArcGap (&helix, &other_helix));
        Hist->deltaT  ->Fill(std::fabs(helix.t0().t0() - other_helix.t0().t0()));
      }
      checked.insert(index); // to avoid double counting in comparisons
    }
    Hist->helicitySum->Fill(helSum);
    Hist->nHelices->Fill(hcol->size());
    Hist->nInSet  ->Fill(indices.size());
  }

  //--------------------------------------------------------------------------------------
  float MultiHelixFilter::evaluateOverlap(const HelixSeed* h1, const HelixSeed* h2) {
    // ensure both exist
    if(!h1 || !h2) return 1.f;

    // get the hit lists for each helix
    std::vector<StrawHitIndex> h1_hits, h2_hits;
    for(size_t index = 0; index < h1->hits().size(); ++index) h1->hits().fillStrawHitIndices(index, h1_hits);
    for(size_t index = 0; index < h2->hits().size(); ++index) h2->hits().fillStrawHitIndices(index, h2_hits);

    // use a set for O(1) searching
    const std::unordered_set<StrawHitIndex> h2_set(h2_hits.begin(), h2_hits.end());

    // check for overlaps
    unsigned overlap = 0;
    for(size_t index = 0; index < h1_hits.size(); ++index) {
      const auto hit = h1_hits[index];
      if(h2_set.find(hit) != h2_set.end()) ++overlap;
    }

    // normalize the overlap to the average number of hits
    const float avg_hits = 0.5 * (h1_hits.size() + h2_hits.size());
    const float ratio = overlap / avg_hits;

    if(_debugLevel > 0) std::cout << ">>> [MultiHelixFilter::" << __func__ << "::" << moduleDescription().moduleLabel() << "]"
                                  << " Overlap = " << overlap << ", Avg hits = " << avg_hits << " ratio = " << ratio
                                  << "\n  Helix 1: " << "N(hits) = " << h1_hits.size() << " p = " << _bz0 * 3./10. * h1->helix().momentum()
                                  << " t = " << h1->t0().t0()
                                  << " x = " << h1->helix().centerx() << " y = " << h1->helix().centery() << " r = " << h1->helix().radius()
                                  << "\n  Helix 2: " << "N(hits) = " << h2_hits.size() << " p = " << _bz0 * 3./10. * h2->helix().momentum()
                                  << " t = " << h2->t0().t0()
                                  << " x = " << h2->helix().centerx() << " y = " << h2->helix().centery() << " r = " << h2->helix().radius()
                                  << std::endl;
    return ratio;
  }

  //--------------------------------------------------------------------------------------
  float MultiHelixFilter::evaluateArcGap(const HelixSeed* h1, const HelixSeed* h2) {
    // ensure both exist
    if(!h1 || !h2) return 0.f;

    // distance between shared circle side
    const float distance_between = std::sqrt(std::pow(h1->helix().centerx() - h2->helix().centerx(), 2) +
                                             std::pow(h1->helix().centery() - h2->helix().centery(), 2));
    const float radial_diff = std::fabs(h1->helix().radius() - h2->helix().radius());
    const float arc_distance = std::fabs(distance_between - radial_diff);

    if(_debugLevel > 0) std::cout << ">>> [MultiHelixFilter::" << __func__ << "::" << moduleDescription().moduleLabel() << "]"
                                  << " distance = " << distance_between << " radial diff = " << radial_diff << " arc_gap = " << arc_distance
                                  << std::endl;
    return arc_distance;
  }

  //--------------------------------------------------------------------------------------
  bool MultiHelixFilter::checkHelixList(std::vector<const HelixSeed*> helices) {
    // check the collection can pass the cuts
    if(helices.empty() || helices.size() > _helixCuts.size()) return false;


    // check if each helix passes the corresponding cut set
    double min_time(helices[0]->t0().t0()), max_time(helices[0]->t0().t0()); // to get the time window of the helices
    for(size_t index = 0; index < helices.size(); ++index) {
      if(!_helixCuts[index].applyCuts(helices[index])) return false;
      min_time = std::min(min_time, helices[index]->t0().t0());
      max_time = std::max(max_time, helices[index]->t0().t0());
    }

    // check the helices are contained in the time window
    if(_timeWindow > 0. && max_time - min_time > _timeWindow) return false;

    // check that all helices pass the overlap check
    const bool do_overlap = _maxOverlap > 0. && _maxOverlap < 1.;
    const bool do_arc     = _minArcGap > 0.;
    if(do_overlap || do_arc) {
      for(size_t i = 0; i < helices.size() - 1; ++i) {
        for(size_t j = i+1; j < helices.size(); ++j) {
          const HelixSeed* h1 = helices[i];
          const HelixSeed* h2 = helices[j];
          if(do_overlap && evaluateOverlap(h1,h2) > _maxOverlap) return false;
          // check if it they share a side of the helix circle
          if(do_arc && evaluateArcGap(h1, h2) < _minArcGap) return false;
        }
      }
    }

    if(_debugLevel > 0) std::cout << ">>> [MultiHelixFilter::" << __func__ << "::" << moduleDescription().moduleLabel() << "]"
                                  << " Accepted a helix set\n";

    // helix set is accepted
    return true;
  }

  //--------------------------------------------------------------------------------------
  void MultiHelixFilter::doPermutations(const std::vector<size_t>& given, const HelixSeedCollection* hcol, std::set<size_t>& accepted) {
    const size_t nhelices = hcol->size();
    if(nhelices < _minNHelices) return; // can't pass the requirements

    if(_debugLevel > 2) {
      std::cout << __func__ << ": Level " << given.size() << ": input set = {";
      for(size_t index = 0; index < given.size(); ++index) {
        if(index > 0) std::cout << ", ";
        std::cout << given[index];
      }
      std::cout << "}, N(helices) = " << nhelices << std::endl;
    }

    // if more than two free indices, continue down
    if(given.size() < _minNHelices) {
      for(size_t index = 0; index < nhelices; ++index) {
        if(find(given.begin(), given.end(), index) != given.end()) continue; //already in the set
        auto next_level = given;
        next_level.push_back(index);
        doPermutations(next_level, hcol, accepted);
      }
      return;
    }

    // at the base level where all indices are fixed
    std::vector<const HelixSeed*> helices;
    for(size_t index = 0; index < given.size(); ++index) helices.push_back(&hcol->at(given[index]));

    // check whether the group passes the cuts
    if(checkHelixList(helices)) {
      for(size_t index = 0; index < given.size(); ++index) {
        const size_t ih = given[index];
        if(accepted.find(ih) == accepted.end()) accepted.insert(ih);
      }
    }
  }

  //--------------------------------------------------------------------------------------
  bool MultiHelixFilter::filter(art::Event& event){
    ++_nevt;

    // create the output
    std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);

    // retrieve the collection
    auto hsH = event.getValidHandle<HelixSeedCollection>(_hsTag);
    const HelixSeedCollection* hscol = hsH.product();

    if(_debugLevel > 1) std::cout << "[MultiHelixFilter::" << __func__ << "::" << moduleDescription().moduleLabel() << "]"
                                  << " Input from: " << _hsTag << " N(helices) = "<< hscol->size() << std::endl;
    if(_debugLevel > 2 && !hscol->empty()) {
      std::cout << "Input helices:\n";
      for(size_t index = 0; index < hscol->size(); ++index) {
        const auto& h = hscol->at(index);
        std::cout << "  Helix " << index << ": " << " p = " << _bz0 * 3./10. * h.helix().momentum()
                  << " t = " << h.t0().t0()
                  << " x = " << h.helix().centerx() << " y = " << h.helix().centery() << " r = " << h.helix().radius()
                  << " E(cluster) = " << ((h.caloCluster().isNull()) ? -1.f : h.caloCluster()->energyDep())
                  << std::endl;
      }
    }

    std::set<size_t> accepted;
    doPermutations({}, hscol, accepted);

    const bool pass = !accepted.empty();

    if(pass) {
      ++_npass;
      if(_debugLevel > 0) std::cout << ">>> [MultiHelixFilter::" << __func__ << "::" << moduleDescription().moduleLabel() << "]"
                                    << " Event " << event.run() << ":" << event.subRun() << ":" << event.event()
                                    << " Accepted an event\n";

      for(auto it = accepted.begin(); it != accepted.end(); ++it) {
        const size_t index = *it;
        if(_debugLevel > 2) std::cout << "  Adding index " << index << " to the output\n";
        triginfo->_helixes.push_back(art::Ptr<HelixSeed>(hsH,index));
      }
    }

    event.put(std::move(triginfo));

    if(_doHistograms) {
      std::set<size_t> all_helices;
      for(size_t index = 0; index < hscol->size(); ++index) all_helices.insert(index);
      fillHistograms(_hists[0], all_helices, hscol);
      if(pass) fillHistograms(_hists[1], accepted, hscol);
    }

    return pass || _noFilter;
  }


  //--------------------------------------------------------------------------------------
  bool MultiHelixFilter::endRun( art::Run& run ) {
    if(_debugLevel > 0) std::cout << "[MultiHelixFilter::" << __func__ << "::" << moduleDescription().moduleLabel() << "]"
                                  << " Passed " <<  _npass << " events out of " << _nevt
                                  << " for a ratio of " << ((_nevt > 0) ? _npass*1./_nevt : 0.)
                                  << std::endl;
    return true;
  }

}
using mu2e::MultiHelixFilter;
DEFINE_ART_MODULE(MultiHelixFilter)
