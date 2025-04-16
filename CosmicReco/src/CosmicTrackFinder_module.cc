// Author : S Middleton
// Data : March 2019
// Purpose: Cosmic Track finder- module calls seed fitting routine to begin cosmic track analysis. The module can call the seed fit and drift fit. Producing a "CosmicTrackSeed" list.

#include "Offline/CosmicReco/inc/CosmicTrackFit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"

//Mu2e General:
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

// ART:
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "art/Utilities/make_tool.h"
#include "canvas/Persistency/Common/Ptr.h"

//MU2E:
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitPosition.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "Offline/TrkReco/inc/TrkTimeCalculator.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/CosmicReco/inc/MinuitDriftFitter.hh"

//utils:
#include "Offline/Mu2eUtilities/inc/ParametricFit.hh"

//For Drift:
#include "Offline/TrkReco/inc/TrkFaceData.hh"

// Mu2e BaBar
#include "Offline/BTrkData/inc/TrkStrawHit.hh"

//CLHEP:
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

//C++:
#include <iostream>
#include <fstream>
#include <string>
#include <memory>
#include <utility>
#include <functional>
#include <float.h>
#include <vector>
#include <map>

using namespace std;
using namespace ROOT::Math::VectorUtil;
using CLHEP::Hep3Vector;
using CLHEP::HepVector;

namespace{
  struct ycomp_iter {
    bool operator()(std::vector<ComboHit>::const_iterator p1, std::vector<ComboHit>::const_iterator p2) { return p1->_pos.y() > p2->_pos.y(); }
  };
}

namespace mu2e {
  class CosmicTrackFinder : public art::EDProducer {
  public:
    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<int>                     debug{Name("debugLevel"), Comment("set to 1 for debug prints"),0};
      fhicl::Atom<int>                     minnsh {Name("minNStrawHits"), Comment("minimum number of straw hits ")};
      fhicl::Atom<int>                     minnch {Name("minNComboHits"), Comment("number of combohits allowed")};
      fhicl::Atom<TrkFitFlag>              saveflag {Name("SaveTrackFlag"),Comment("if set to OK then save the track"),TrkFitFlag::helixOK};
      fhicl::Atom<int>                     minNHitsTimeCluster{Name("minNHitsTimeCluster"),Comment("minium allowed time cluster")};
      fhicl::Atom<art::InputTag>           chToken{Name("ComboHitCollection"),Comment("tag for combo hit collection")};
      fhicl::Atom<art::InputTag>           tcToken{Name("TimeClusterCollection"),Comment("tag for time cluster collection")};
      fhicl::Atom<bool>                    UseLineFinder{Name("UseLineFinder"),Comment("use line finder for seeding drift fit")};
      fhicl::Atom<bool>                    UseChiFit{Name("UseChiFit"),Comment("use chi fit to improve seed")};
      fhicl::Atom<art::InputTag>           lfToken{Name("LineFinderTag"),Comment("tag for line finder seed")};
      fhicl::Atom<bool>                    DoDrift{Name("DoDrift"),Comment("turn on for drift fit")};
      fhicl::Atom<bool>                    UseTime{Name("UseTime"),Comment("use time for drift fit")};
      fhicl::Atom<double>                  driftRes{Name("DriftRes"),Comment("Drift resolution for first fit stage")};
      fhicl::Atom<double>                  mnTolerance{Name("MinuitTolerance"),Comment("Tolerance for minuit convergence")};
      fhicl::Atom<double>                  mnPrecision{Name("MinuitPrecision"),Comment("Effective precision for likelihood function")};
      fhicl::Table<CosmicTrackFit::Config> tfit{Name("CosmicTrackFit"), Comment("fit")};
    };
    typedef art::EDProducer::Table<Config> Parameters;

    explicit CosmicTrackFinder(const Parameters& conf);
    virtual ~CosmicTrackFinder();

    virtual void beginJob() override;
    virtual void beginRun(art::Run& run) override;
    virtual void produce(art::Event& event ) override;

  private:

    Config _conf;

    int                                 _debug;
    int                                 _minnsh; // minimum # of strawHits in CH
    int                                 _minnch; // minimum # of ComboHits for viable fit
    TrkFitFlag        _saveflag;//write tracks that satisfy these flags
    int                                 _minNHitsTimeCluster; //min number of hits in a time cluster
    //float                                _max_seed_chi2; ///maximum chi2 allowed for seed

    art::InputTag  _chToken;
    art::InputTag  _tcToken;

    bool _UseLineFinder;
    bool _UseChiFit;
    art::InputTag _lfToken;

    bool            _DoDrift;
    bool       _UseTime;
    double _driftRes;
    double _mnTolerance;
    double _mnPrecision;

    CosmicTrackFit     _tfit;

    ProditionsHandle<StrawResponse> _strawResponse_h;
    ProditionsHandle<Tracker> _alignedTracker_h;

    void     OrderHitsY(ComboHitCollection const&chcol, std::vector<StrawHitIndex> const&inputIdx, std::vector<StrawHitIndex> &outputIdxs);
    int      goodHitsTimeCluster(const TimeCluster &TCluster, ComboHitCollection const& chcol);

  };
  CosmicTrackFinder::CosmicTrackFinder(const Parameters& conf) :
    art::EDProducer(conf),
    _debug  (conf().debug()),
    _minnsh   (conf().minnsh()),
    _minnch  (conf().minnch()),
    _saveflag  (conf().saveflag()),
    _minNHitsTimeCluster(conf().minNHitsTimeCluster()),
    _chToken (conf().chToken()),
    _tcToken (conf().tcToken()),
    _UseLineFinder (conf().UseLineFinder()),
    _UseChiFit (conf().UseChiFit()),
    _lfToken (conf().lfToken()),
    _DoDrift (conf().DoDrift()),
    _UseTime (conf().UseTime()),
    _driftRes(conf().driftRes()),
    _mnTolerance (conf().mnTolerance()),
    _mnPrecision (conf().mnPrecision()),
    _tfit (conf().tfit())
  {
    consumes<ComboHitCollection>(_chToken);
    consumes<TimeClusterCollection>(_tcToken);
    mayConsume<CosmicTrackSeedCollection>(_lfToken);
    produces<CosmicTrackSeedCollection>();

  }

  CosmicTrackFinder::~CosmicTrackFinder(){}

  void CosmicTrackFinder::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
  }

  void CosmicTrackFinder::beginRun(art::Run& run) {
  }

  void CosmicTrackFinder::produce(art::Event& event ) {
    Tracker const& tracker = _alignedTracker_h.get(event.id());
    _tfit.setTracker(&tracker);
    StrawResponse const& srep = _strawResponse_h.get(event.id());

    if (_debug != 0) {
      std::cout << "CosmicTrackFinder: Producing Cosmic Track ..."<<std::endl;
    }

    unique_ptr<CosmicTrackSeedCollection> seed_col = make_unique<CosmicTrackSeedCollection>();

    int _iev = event.id().event();

    if (_debug > 0){
      std::cout << "CosmicTrackFinder: ST Finder Event #" << _iev << std::endl;
    }

    auto const& chH = event.getValidHandle<ComboHitCollection>(_chToken);
    const ComboHitCollection& chcol(*chH);
    auto  const& tcH = event.getValidHandle<TimeClusterCollection>(_tcToken);
    const TimeClusterCollection& tccol(*tcH);

    for (size_t index=0;index< tccol.size();++index) {
      int   nGoodTClusterHits(0);
      const auto& tclust = tccol[index];
      nGoodTClusterHits     = goodHitsTimeCluster(tclust,chcol);

      if ( nGoodTClusterHits < _minNHitsTimeCluster) {  continue; }
      if (_debug > 0){
        std::cout<<"CosmicTrackFinder: time clusters " << _iev << std::endl;
      }

      std::vector<StrawHitIndex> panelHitIdxs;
      OrderHitsY(chcol,tclust.hits(),panelHitIdxs);

      int nFiltComboHits = 0;
      int nFiltStrawHits = 0;
      for (size_t i=0;i<panelHitIdxs.size();i++){
        auto ch = chcol[panelHitIdxs[i]];
        nFiltComboHits++;
        nFiltStrawHits += ch.nStrawHits();
      }

      if (_debug != 0){
        std::cout<<"CosmicTrackFinder: #filtered SHits"<<nFiltStrawHits<<" #filter CHits "<<nFiltComboHits<<std::endl;
      }

      if (nFiltComboHits < _minnch ) {        continue; }
      if (nFiltStrawHits < _minnsh) { continue; }


      ostringstream title;
      title << "Run: " << event.id().run()
            << "  Subrun: " << event.id().subRun()
            << "  Event: " << event.id().event()<<".root";

      CosmicTrackSeed tseed ;
      if (_UseLineFinder) {
        auto const& lfH = event.getValidHandle<CosmicTrackSeedCollection>(_lfToken);
        const CosmicTrackSeedCollection& lfcol(*lfH);
        if (lfcol.size() == 0) {
          continue;
        }

        tseed = lfcol[0];
        double _interror = 40;
        double _direrror = 2.5;
        double _t0error = 1;

        tseed._track.FitParams.Covarience.sigA0 = _interror;
        tseed._track.FitParams.Covarience.sigA1 = _direrror;
        tseed._track.FitParams.Covarience.sigB0 = _interror;
        tseed._track.FitParams.Covarience.sigB1 = _direrror;
        tseed._t0._t0err = _t0error;

      } else{
        tseed._t0          = tclust._t0;
        tseed._timeCluster = art::Ptr<TimeCluster>(tcH,index);
        tseed._status.merge(TrkFitFlag::Straight);
        tseed._status.merge(TrkFitFlag::hitsOK);
      }

      if (_UseChiFit){
        _tfit.BeginFit(title.str().c_str(), tseed, event, chcol, panelHitIdxs);

        if( _tfit.goodTrack(tseed._track) == false){
          tseed._status.clear(TrkFitFlag::helixConverged);
          tseed._status.clear(TrkFitFlag::helixOK);
        }
      }

      if (tseed._status.hasAnyProperty(TrkFitFlag::helixOK) && tseed._status.hasAnyProperty(TrkFitFlag::helixConverged) && tseed._track.converged == true ) {

        if (tseed.status().hasAnyProperty(_saveflag)){

          if(_DoDrift) {
            if (_UseTime) {
              MinuitDriftFitter::DoDriftTimeFit(_debug,tseed, srep, &tracker, _driftRes, _mnTolerance, _mnPrecision );
            } else {
              _tfit.DriftFit(tseed, srep);
            }

            if( !tseed._track.minuit_converged ){
              continue;
            }
            // TODO: do we keep or remove this?
            // // remove all outliers from ComboHitCollection
            // ComboHitCollection tmpHits;
            // for(auto const &chit : tseed._straw_chits){
            //   if(!chit._flag.hasAnyProperty(StrawHitFlag::outlier)){
            //     tmpHits.push_back(chit);
            //   }
            // }
            // tseed._straw_chits = tmpHits;
            // if (tmpHits.size() == 0)
            //   continue;
          }else{
            tseed._track.MinuitParams = tseed._track.FitParams;
            tseed._track.MinuitCoordSystem = tseed._track.FitCoordSystem;
            tseed._track.MinuitEquation = tseed._track.FitEquation;
            tseed._track.MinuitParams.cov = std::vector<double>(15, 0);
          }

          tseed._straw_chits.setParent(chcol.parent());

          CosmicTrackSeedCollection* col = seed_col.get();
          col->emplace_back(tseed);
        }
      }
    }

    event.put(std::move(seed_col));
  }

  void CosmicTrackFinder::OrderHitsY(ComboHitCollection const& chcol, std::vector<StrawHitIndex> const& inputIdxs, std::vector<StrawHitIndex> &outputIdxs){
    if (_debug != 0){
      std::cout<<"Ordering Hits..."<<std::endl;
    }

    std::vector<std::vector<ComboHit>::const_iterator> ordChColIters;
    for (size_t i=0;i<inputIdxs.size();i++){
      std::vector<ComboHit>::const_iterator tempiter = chcol.begin() + inputIdxs[i];
      ordChColIters.push_back(tempiter);
    }
    std::sort(ordChColIters.begin(), ordChColIters.end(),ycomp_iter());

    for (size_t i=0;i<ordChColIters.size();i++){
      auto thisiter = ordChColIters[i];
      outputIdxs.push_back(std::distance(chcol.begin(),thisiter));
    }
  }

  int  CosmicTrackFinder::goodHitsTimeCluster(TimeCluster const& TCluster, ComboHitCollection const& chcol){
    int nhits = TCluster.nhits();
    int ngoodhits(0);
    double minT(500.), maxT(2000.);

    for (int i=0; i<nhits; ++i) {
      int          index   = TCluster.hits().at(i);
      ComboHit     sh      = chcol.at(index);
      if ( (sh.time() < minT) || (sh.time() > maxT) )  { continue; }

      ngoodhits += sh.nStrawHits();
    }

    return ngoodhits;
  }

}
using mu2e::CosmicTrackFinder;
DEFINE_ART_MODULE(CosmicTrackFinder)
