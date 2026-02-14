//
// A Filter module aimed to select events using a Likelihood defined with calorimeter cluster info
//
//
// Original author G. Pezzullo
//

#include "CLHEP/Units/SystemOfUnits.h"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"

#include "Offline/ConfigTools/inc/ConfigFileLookupPolicy.hh"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/CalorimeterGeom/inc/DiskCalorimeter.hh"

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"

#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"

// #include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Provenance.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/Selector.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"

#include "canvas/Utilities/InputTag.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"
#include "fhiclcpp/types/Sequence.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include <cmath>
#include <string>
#include <vector>

namespace mu2e {

class CaloClusterCounter : public art::EDFilter {

public:
  using Name = fhicl::Name;
  using Comment = fhicl::Comment;
  struct Config {
    fhicl::Atom<int> diag{Name("diagLevel"), Comment("Diagnostic Level"), 0};
    fhicl::Atom<art::InputTag> CCTag{Name("CaloClusterModuleLabel"),
                                     Comment("CaloClusterModuleLabel producer")};
    fhicl::Atom<double> minClE{Name("MinClusterEnergy"),
                               Comment("Minimum energy deposit in a calo cluster")};
    fhicl::Atom<double> maxClE{Name("MaxClusterEnergy"),
                               Comment("Maximum energy deposit in a calo cluster")};
    fhicl::Atom<double> minClR{Name("MinClusterRadius"),
                               Comment("Minimum radial distance of the cluster from the DS axis")};
    fhicl::Atom<int> minNCl{Name("MinNCl"), Comment("Minimum number of calo cluster")};
    fhicl::Sequence<int> diskID{Name("DiskID"), Comment("Disk ID")};
    fhicl::Atom<int> minNCel{Name("MinNCel"),
                             Comment("Minimum number of crystals within the cluster")};
    fhicl::Atom<int> maxNCel{Name("MaxNCel"),
                             Comment("Maximum number of crystals within the cluster")};
  };

  virtual ~CaloClusterCounter() {}

  virtual void beginJob();
  virtual void endJob();
  virtual bool filter(art::Event& event) override;
  virtual bool beginRun(art::Run& run);
  virtual bool endRun(art::Run& run) override;

  using Parameters = art::EDFilter::Table<Config>;
  explicit CaloClusterCounter(const Parameters& conf);

private:
  int _diagLevel;
  int _nProcess;
  int _nPass;
  art::InputTag _clTag;
  double _minClEnergy;
  double _maxClEnergy;
  double _minClRadius;
  int _minNCl;
  std::vector<int> _diskID;
  int _minNCel;
  int _maxNCel;

  const Calorimeter* _calogeom;
};

CaloClusterCounter::CaloClusterCounter(const Parameters& config) :
    art::EDFilter{config}, _diagLevel(config().diag()), _nProcess(0), _nPass(0),
    _clTag(config().CCTag()), _minClEnergy(config().minClE()), _maxClEnergy(config().maxClE()),
    _minClRadius(config().minClR()), _minNCl(config().minNCl()), _diskID(config().diskID()),
    _minNCel(config().minNCel()), _maxNCel(config().maxNCel()) {

  produces<TriggerInfo>();
}

void CaloClusterCounter::beginJob() {}

void CaloClusterCounter::endJob() {}

bool CaloClusterCounter::beginRun(art::Run& run) {
  GeomHandle<Calorimeter> ch;
  _calogeom = ch.get();

  return true;
}

bool CaloClusterCounter::endRun(art::Run& run) {
  if (_diagLevel > 0 && _nProcess > 0) {
    std::cout << "[CaloClusterCounter::endRun]"
              << " passed " << _nPass << " events out of " << _nProcess << " for a ratio of "
              << float(_nPass) / float(_nProcess) << std::endl;
  }
  return true;
}

//--------------------------------------------------------------------------------
// Follow the body of the Filter logic
//--------------------------------------------------------------------------------
bool CaloClusterCounter::filter(art::Event& event) {

  ++_nProcess;
  if (_nProcess % 10 == 0 && _diagLevel > 0)
    std::cout << "[CaloClusterCounter::filter] Processing event from CaloClusterCounter =  "
              << _nProcess << std::endl;

  std::unique_ptr<TriggerInfo> triginfo(new TriggerInfo);
  bool retval(false);

  // Get calo cluster collection
  auto clH = event.getValidHandle<CaloClusterCollection>(_clTag);
  const CaloClusterCollection* caloClusters = clH.product();

  int nClusterAboveThreshold(0);

  // for loop over the clusters in the calorimeter
  for (auto icl = caloClusters->begin(); icl != caloClusters->end(); ++icl) {
    auto const& cluster = *icl;
    CLHEP::Hep3Vector pos = cluster.cog3Vector();
    float clRadius = std::sqrt(std::pow(pos.x(), 2) + std::pow(pos.y(), 2));
    if (_diagLevel > 0) {
      std::cout << "[CaloClusterCounter::filter] cluster "
                << std::distance(caloClusters->begin(), icl) << " diskID = " << cluster.diskID()
                << " E = " << cluster.energyDep() << " radial distance = " << clRadius
                << " nCels = " << cluster.size() << std::endl;
    }
    if ((std::find(_diskID.begin(), _diskID.end(), cluster.diskID()) != _diskID.end()) &&
        (cluster.energyDep() >= _minClEnergy) && (cluster.energyDep() <= _maxClEnergy) &&
        (cluster.size() >= _minNCel) && (cluster.size() <= _maxNCel) &&
        (clRadius >= _minClRadius)) {
      size_t index = std::distance(caloClusters->begin(), icl);
      triginfo->_caloClusters.push_back(art::Ptr<CaloCluster>(clH, index));
      if(_diagLevel > 0) std::cout << "  --> Accepted cluster\n";

      ++nClusterAboveThreshold;
    } else if(_diagLevel > 1) std::cout << "  --> Rejected cluster\n";

  }

  if (nClusterAboveThreshold >= _minNCl) {
    if(_diagLevel > 0) std::cout << "[CaloClusterCounter::filter] Event " << event.id()
                                 << ": Accepting event, N(accepted clusters) = " << nClusterAboveThreshold
                                 << " from " << caloClusters->size() << " clusters"
                                 << std::endl;
    retval = true;
  } else if(_diagLevel > 1) std::cout << "[CaloClusterCounter::filter] Event " << event.id()
                                      << ": Rejecting event, N(accepted clusters) = " << nClusterAboveThreshold
                                      << " from " << caloClusters->size() << " clusters"
                                      << std::endl;


  event.put(std::move(triginfo));
  return retval;
}

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::CaloClusterCounter)
