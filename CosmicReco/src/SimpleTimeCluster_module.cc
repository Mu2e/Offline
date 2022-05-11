//Author: R Bonventre
//Date: Feb 2020
//Purpose: To improve TimeClustering for Cosmics
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/EDProducer.h"
#include "art_root_io/TFileService.h"
// Mu2e
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"
// data
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
// tracking
#include "Offline/TrkReco/inc/TrkUtilities.hh"
#include "Offline/TrkReco/inc/TrkTimeCalculator.hh"
// C++
#include <memory>
#include <algorithm>
#include <utility>
using namespace std;

namespace {
  // comparison functor for sorting by t
  struct tcomp : public std::binary_function<mu2e::ComboHit,mu2e::ComboHit,bool> {
    bool operator()(mu2e::ComboHit const& p1, mu2e::ComboHit const& p2) { return p1.time() < p2.time(); }
  };
}



namespace mu2e {
  class SimpleTimeCluster : public art::EDProducer {
    public:
      struct Config{
        using Name=fhicl::Name;
        using Comment=fhicl::Comment;
        fhicl::Atom<int> debug{Name("debugLevel"), Comment("set to 1 for debug prints"),0};
        fhicl::Atom<int> minnsh {Name("minNStrawHits"), Comment("minimum number of straw hits "),5};
        fhicl::OptionalAtom<int> maxnsh {Name("maxNStrawHits"), Comment("maximum number of straw hits ")};
        fhicl::Atom<double> timewindow {Name("TimeWindow"), Comment("Width of time window in ns"),100};
        fhicl::Atom<bool> testflag {Name("TestFlag"),Comment("Test StrawHitFlags")};
        fhicl::Atom<bool> useonepanel {Name("UseOnlyOnePanel"),Comment("Use hits from one panel only"),false};
        fhicl::Atom<art::InputTag> chToken{Name("ComboHitCollection"),Comment("tag for combo hit collection")};
        fhicl::Atom<art::InputTag> shfToken{Name("StrawHitFlagCollection"),Comment("tag for StrawHitFlag collection")};
        fhicl::Sequence<std::string> hsel{Name("HitSelectionBits"),Comment("Required flags if TestFlag is set"),vector<string>{"EnergySelection","TimeSelection"}};
        fhicl::Sequence<std::string> hbkg{Name("HitBackgroundBits"),Comment("Excluded flags if TestFlag is set"),vector<string>{}};
      };
      typedef art::EDProducer::Table<Config> Parameters;
      explicit SimpleTimeCluster(const Parameters& conf);

      void produce(art::Event& e) override;

    private:
      int               _iev;
      int               _debug;
      int               _minnsh;
      bool              _hasmaxnsh;
      int               _maxnsh;
      double            _timeWindow;
      bool              _testflag;
      bool _useonepanel;
      art::InputTag  _chToken;
      art::InputTag _shfToken;
      const StrawHitFlagCollection *_shfcol;
      const ComboHitCollection *_chcol;
      StrawHitFlag      _hsel, _hbkg;


      void findClusters(TimeClusterCollection& tccol);
      bool goodHit(const StrawHitFlag& flag) const;
  };

  SimpleTimeCluster::SimpleTimeCluster(const Parameters& conf) :
    art::EDProducer(conf),
    _debug (conf().debug()),
    _minnsh (conf().minnsh()),
    _hasmaxnsh (false),
    _maxnsh (0),
    _timeWindow (conf().timewindow()),
    _testflag (conf().testflag()),
    _useonepanel (conf().useonepanel()),
    _chToken (conf().chToken()),
    _shfToken (conf().shfToken()),
    _hsel (conf().hsel()),
    _hbkg (conf().hbkg())
  {
    _hasmaxnsh = conf().maxnsh(_maxnsh);
    produces<TimeClusterCollection>();
  }


  //--------------------------------------------------------------------------------------------------------------
  void SimpleTimeCluster::produce(art::Event & event ){
    _iev = event.id().event();

    auto const& chH = event.getValidHandle<ComboHitCollection>(_chToken);
    _chcol = chH.product();

    if(_testflag){
      auto shfH = event.getValidHandle<StrawHitFlagCollection>(_shfToken);
      _shfcol = shfH.product();
      if(_shfcol->size() != _chcol->size())
        throw cet::exception("RECO")<<"SimpleTimeCluster: inconsistent flag collection length " << endl;
    }

    std::unique_ptr<TimeClusterCollection> tccol(new TimeClusterCollection);
    findClusters(*tccol);

    if (_debug > 0) std::cout << "Found " << tccol->size() << " Time Clusters " << std::endl;

    if (_debug > 1){
      for(auto const& tc : *tccol) {
        std::cout << "Time Cluster time = " << tc.t0().t0() << " +- " << tc.t0().t0Err()
          << " position = " << tc._pos << std::endl;
        if(_debug > 3){
          for (auto shi : tc._strawHitIdxs ) {
            std::cout << "Time Cluster hit at index " << shi << std::endl;
          }
        }
      }
    }
    event.put(std::move(tccol));
  }

  void SimpleTimeCluster::findClusters(TimeClusterCollection& tccol) {
    //sort the hits by time
    ComboHitCollection ordChCol;
    ordChCol.reserve(_chcol->size());

    for (size_t i=0; i<_chcol->size(); ++i) {
      const ComboHit& ch  = _chcol->at(i);
      if (_testflag && !goodHit((*_shfcol)[i]))
        continue;
      ordChCol.push_back(ComboHit(ch));
    }
    if (ordChCol.size() == 0) return;

    std::sort(ordChCol.begin(), ordChCol.end(),tcomp());

    size_t maxStart = 0;
    size_t maxEnd = 0;
    int maxCount = 1;

    size_t endIndex = 1;
    int count = ordChCol[0].nStrawHits();
    for (size_t startIndex = 0;startIndex < ordChCol.size();startIndex++){
      double startTime = ordChCol[startIndex].time();
      while (true){
        if (endIndex >= ordChCol.size())
          break;
        if (ordChCol[endIndex].time()-startTime > _timeWindow)
          break;
        count += ordChCol[endIndex].nStrawHits();
        if (count > maxCount){
          maxCount = count;
          maxStart = startIndex;
          maxEnd = endIndex;
        }
        endIndex++;
      }
      count -= ordChCol[startIndex].nStrawHits();
    }

    if (maxCount < _minnsh) return;
    if (_hasmaxnsh && maxCount > _maxnsh) return;

    std::vector<int> hits_in_panel(StrawId::_nupanels,0);
    if (_useonepanel){
      for (size_t i=0; i<_chcol->size(); ++i) {
        if (_testflag && !goodHit((*_shfcol)[i]))
          continue;
        if (_chcol->at(i).time() < ordChCol[maxStart].time() || _chcol->at(i).time() > ordChCol[maxEnd].time())
          continue;
        hits_in_panel[_chcol->at(i).strawId().uniquePanel()] += _chcol->at(i).nStrawHits();
      }
    }
    size_t max_panel = std::distance(hits_in_panel.begin(),std::max_element(hits_in_panel.begin(),hits_in_panel.end()));


    TimeCluster tclust;
    tclust._nsh = 0;
    double avg = 0;
    for (size_t i=0;i<_chcol->size();i++){
      if (_testflag && !goodHit((*_shfcol)[i]))
        continue;
      if (_chcol->at(i).time() < ordChCol[maxStart].time() || _chcol->at(i).time() > ordChCol[maxEnd].time())
        continue;
      if (_useonepanel && (_chcol->at(i).strawId().uniquePanel() != max_panel))
        continue;
      avg += _chcol->at(i).time();
      tclust._strawHitIdxs.push_back(i);
      tclust._nsh += _chcol->at(i).nStrawHits();
    }
    tclust._t0 = TrkT0(avg/tclust._strawHitIdxs.size(),(ordChCol[maxEnd].time()-ordChCol[maxStart].time())/2.);
    tccol.push_back(tclust);
  }

  bool SimpleTimeCluster::goodHit(const StrawHitFlag& flag) const
  {
    return flag.hasAllProperties(_hsel) && !flag.hasAnyProperty(_hbkg);
  }

}

using mu2e::SimpleTimeCluster;
DEFINE_ART_MODULE(SimpleTimeCluster);
