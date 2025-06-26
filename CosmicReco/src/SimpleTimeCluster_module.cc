// Author: R Bonventre
// Date: Feb 2020
// Purpose: To improve TimeClustering for Cosmics
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
// Mu2e
#include "Offline/GeneralUtilities/inc/Angles.hh"
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "Offline/Mu2eUtilities/inc/polyAtan2.hh"
// data
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
// tracking
#include "Offline/TrkReco/inc/TrkTimeCalculator.hh"
#include "Offline/TrkReco/inc/TrkUtilities.hh"
// C++
#include <algorithm>
#include <memory>
#include <utility>
#include <numeric>

namespace {
// comparison functor for sorting by t
struct tcomp {
  bool operator()(mu2e::ComboHit const& p1, mu2e::ComboHit const& p2) {
    return p1.correctedTime() < p2.correctedTime();
  }
};
} // namespace

namespace mu2e {
class SimpleTimeCluster : public art::EDProducer {
public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;
    fhicl::Atom<int> debug{Name("debugLevel"), Comment("set to 1 for debug prints")};
    fhicl::Atom<int> minnsh{Name("minNStrawHits"), Comment("minimum number of straw hits ")};
    fhicl::Atom<int> minnpanels{Name("minNPanels"), Comment("minimum number of panels ")};
    fhicl::OptionalAtom<int> maxnsh{Name("maxNStrawHits"), Comment("maximum number of straw hits ")};
    fhicl::Atom<int> timewindow{Name("TimeWindow"), Comment("Width of time window in ns")};
    fhicl::Atom<bool> usetimewindow{Name("useTimeWindow"), Comment("Use timewindow cut")};
    fhicl::Atom<int> timestep{Name("TimeStep"), Comment("Width of timestep between hits in ns")};
    fhicl::Atom<bool> usetimestep{Name("useTimeStep"), Comment("Use timestep cut"), true};
    fhicl::Atom<bool> testflag{Name("TestFlag"), Comment("Test StrawHitFlags")};
    fhicl::Atom<bool> useonepanel{Name("UseOnlyOnePanel"), Comment("Use hits from one panel only")};
    fhicl::Atom<bool> useoneplane{Name("UseOnlyOnePlane"), Comment("Use hits from one plane only")};
    fhicl::Atom<art::InputTag> chToken{Name("ComboHitCollection"), Comment("tag for combo hit collection")};
    fhicl::Sequence<std::string> hsel{Name("HitSelectionBits"), Comment("Required flags if TestFlag is set")};
    fhicl::Sequence<std::string> hbkg{Name("HitBackgroundBits"), Comment("Excluded flags if TestFlag is set")};
    fhicl::Sequence<std::string> hnotnoise{Name("HitNonNoiseSelectionBits"), Comment("Require at least one hit with these flags")};
  };
  typedef art::EDProducer::Table<Config> Parameters;
  explicit SimpleTimeCluster(const Parameters& conf);

  void produce(art::Event& e) override;

private:
  int _iev;
  int _debug;
  int _minnsh;
  int _minnpanels;
  bool _hasmaxnsh;
  int _maxnsh;
  bool _usetimeWindow;
  int _timeWindow;
  bool _usetimeStep;
  int _timeStep;
  bool _testflag;
  bool _useonepanel;
  bool _useoneplane;
  art::InputTag _chToken;
  const ComboHitCollection* _chcol;
  StrawHitFlag _hsel, _hbkg, _hnotnoise;

  void findClusters(TimeClusterCollection& tccol);
  bool goodHit(const StrawHitFlag& flag) const;
};

SimpleTimeCluster::SimpleTimeCluster(const Parameters& conf) :
    art::EDProducer(conf),
    _debug(conf().debug()),
    _minnsh(conf().minnsh()),
    _minnpanels(conf().minnpanels()),
    _hasmaxnsh(false),
    _maxnsh(0),
    _usetimeWindow(conf().usetimewindow()),
    _timeWindow(conf().timewindow()),
    _usetimeStep(conf().usetimestep()),
    _timeStep(conf().timestep()),
    _testflag(conf().testflag()),
    _useonepanel(conf().useonepanel()),
    _useoneplane(conf().useoneplane()),
    _chToken(conf().chToken()),
    _hsel(conf().hsel()),
    _hbkg(conf().hbkg()),
    _hnotnoise(conf().hnotnoise())
{
  _hasmaxnsh = conf().maxnsh(_maxnsh);
  produces<TimeClusterCollection>();
}

//--------------------------------------------------------------------------------------------------------------
void SimpleTimeCluster::produce(art::Event& event) {
  _iev = event.id().event();

  auto const& chH = event.getValidHandle<ComboHitCollection>(_chToken);
  _chcol = chH.product();

  std::unique_ptr<TimeClusterCollection> tccol(new TimeClusterCollection);
  findClusters(*tccol);

  if (_debug > 0)
    std::cout << "Found " << tccol->size() << " Time Clusters " << std::endl;

  if (_debug > 1) {
    for (auto const& tc : *tccol) {
      std::cout << "Time Cluster time = " << tc.t0().t0() << " +- " << tc.t0().t0Err()
                << " position = " << tc._pos << " nhits = " << tc.hits().size() << std::endl;
      if (_debug > 3) {
        for (auto shi : tc._strawHitIdxs) {
          std::cout << "Time Cluster hit at index " << shi << std::endl;
        }
      }
    }
  }
  event.put(std::move(tccol));
}

void SimpleTimeCluster::findClusters(TimeClusterCollection& tccol) {
  // sort the hits by time
  ComboHitCollection ordChCol;
  ordChCol.reserve(_chcol->size());

  for (size_t i = 0; i < _chcol->size(); ++i) {
    const ComboHit& ch = _chcol->at(i);
    if (_testflag && !goodHit(ch.flag()))
      continue;
    ordChCol.push_back(ComboHit(ch));
  }
  if (ordChCol.size() == 0)
    return;

  std::sort(ordChCol.begin(), ordChCol.end(), tcomp());

  // modified to search for multiple peaks
  std::vector<size_t> peakStart, peakEnd;
  std::vector<size_t> peakHits;
  size_t endIndex = 0;
  int count = 0;
  for (size_t startIndex = 0; startIndex < ordChCol.size(); startIndex++) {
    count = ordChCol[startIndex].nStrawHits();
    bool hasNonNoise = false;
    if (endIndex >= ordChCol.size())
      break;
    if (startIndex < endIndex)
      continue;
    if (startIndex > endIndex)
      endIndex = startIndex;
    double startTime = ordChCol[startIndex].correctedTime();
    double endTime = -1;
    if (ordChCol[endIndex].flag().hasAnyProperty(_hnotnoise))
      hasNonNoise = true;
    while (true) {
      endIndex++;
      if (endIndex >= ordChCol.size())
        break;
      endTime = ordChCol[endIndex].correctedTime();
      count += ordChCol[endIndex].nStrawHits();
      if (ordChCol[endIndex].flag().hasAnyProperty(_hnotnoise))
        hasNonNoise = true;
      if (_usetimeStep && endTime - ordChCol[endIndex - 1].correctedTime() > _timeStep)
        break;
      if (endTime - startTime < 0)
        break;
      if (_usetimeWindow && endTime - startTime > _timeWindow)
        break;
    }
    if (endIndex < ordChCol.size())
      count -= ordChCol[endIndex].nStrawHits();
    endIndex--;
    endTime = ordChCol[endIndex].correctedTime();
    if (_hasmaxnsh && count > _maxnsh)
      continue;
    if (!_usetimeWindow && endTime - startTime > _timeWindow)
      continue;
    if (count >= _minnsh && hasNonNoise) {
      peakStart.push_back(startIndex);
      peakEnd.push_back(endIndex);
      peakHits.push_back(count);
    }
  }
  if (peakHits.size() == 0)
    return;

  std::vector<size_t> max_panels(peakHits.size(),-1);
  std::vector<size_t> max_planes(peakHits.size(),-1);
  for (size_t n = 0; n < peakHits.size(); n++) {
    std::vector<int> hits_in_panel(StrawId::_nupanels, 0);
    std::vector<int> hits_in_plane(StrawId::_nplanes, 0);
    if (_useonepanel) {
      for (size_t i = 0; i < _chcol->size(); ++i) {
        if (_testflag && !goodHit(_chcol->at(i).flag()))
          continue;
        if (_chcol->at(i).correctedTime() < ordChCol[peakStart[n]].correctedTime() ||
            _chcol->at(i).correctedTime() > ordChCol[peakEnd[n]].correctedTime())
          continue;
        hits_in_panel[_chcol->at(i).strawId().uniquePanel()] += _chcol->at(i).nStrawHits();
        hits_in_plane[_chcol->at(i).strawId().plane()] += _chcol->at(i).nStrawHits();
      }
    }

    int npanels = 0;
    for (size_t ip = 0; ip < hits_in_panel.size(); ip++) {
      if (hits_in_panel[ip] > 0)
        npanels++;
    }
    if (_minnpanels > 0 && npanels < _minnpanels)
      continue;
    max_panels[n] = std::distance(hits_in_panel.begin(), std::max_element(hits_in_panel.begin(), hits_in_panel.end()));
    max_planes[n] = std::distance(hits_in_plane.begin(), std::max_element(hits_in_plane.begin(), hits_in_plane.end()));
  }

  // Record TimeCluster info
  size_t nClust = peakHits.size();
  double avg = 0;
  int chCount;
  for (size_t n = 0; n < nClust; n++) {
    std::vector<int> hits_in_panel(StrawId::_nupanels, 0);
    TimeCluster tclust;
    tclust._nsh = peakHits[n];
    avg = 0;
    chCount = 0;
    double time1 = ordChCol[peakStart[n]].correctedTime();
    double time2 = ordChCol[peakEnd[n]].correctedTime();
    for (size_t i = 0; i < _chcol->size(); i++) {

      if (_testflag && !goodHit(_chcol->at(i).flag()))
        continue;
      if (_chcol->at(i).correctedTime() < time1 || _chcol->at(i).correctedTime() > time2)
        continue;
      if (_useonepanel && (_chcol->at(i).strawId().uniquePanel() != max_panels[n]))
        continue;
      if (_useoneplane && (_chcol->at(i).strawId().plane() != max_planes[n]))
        continue;
      hits_in_panel[_chcol->at(i).strawId().uniquePanel()] = 1;
      avg += _chcol->at(i).correctedTime();
      tclust._strawHitIdxs.push_back(i);
      chCount++;
    }
    int npanels = std::accumulate(hits_in_panel.begin(),hits_in_panel.end(),0);
    if (_minnpanels > 0 && npanels < _minnpanels)
      continue;

    tclust._t0 = TrkT0(avg / chCount, (ordChCol[peakEnd[n]].correctedTime() - ordChCol[peakStart[n]].correctedTime()) / 2.);
    tccol.push_back(tclust);
  }
}

bool SimpleTimeCluster::goodHit(const StrawHitFlag& flag) const {
  return flag.hasAllProperties(_hsel) && !flag.hasAnyProperty(_hbkg);
}

} // namespace mu2e

using mu2e::SimpleTimeCluster;
DEFINE_ART_MODULE(SimpleTimeCluster)
