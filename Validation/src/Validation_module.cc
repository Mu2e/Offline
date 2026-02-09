//
// This module loops over the data products in the event
// and makes a set of standard validation histograms
// for instances of products it knows about
//
// Ray Culbertson
//

#include "Offline/Validation/inc/ValBkgCluster.hh"
#include "Offline/Validation/inc/ValBkgQual.hh"
#include "Offline/Validation/inc/ValCaloCluster.hh"
#include "Offline/Validation/inc/ValCaloDigi.hh"
#include "Offline/Validation/inc/ValCaloHit.hh"
#include "Offline/Validation/inc/ValCaloRecoDigi.hh"
#include "Offline/Validation/inc/ValCaloShowerStep.hh"
#include "Offline/Validation/inc/ValComboHit.hh"
#include "Offline/Validation/inc/ValCrvCoincidenceCluster.hh"
#include "Offline/Validation/inc/ValCrvDigi.hh"
#include "Offline/Validation/inc/ValCrvDigiMC.hh"
#include "Offline/Validation/inc/ValCrvRecoPulse.hh"
#include "Offline/Validation/inc/ValCrvStep.hh"
#include "Offline/Validation/inc/ValEventWindowMarker.hh"
#include "Offline/Validation/inc/ValGenParticle.hh"
#include "Offline/Validation/inc/ValHelixSeed.hh"
#include "Offline/Validation/inc/ValKalSeed.hh"
#include "Offline/Validation/inc/ValProtonBunchIntensity.hh"
#include "Offline/Validation/inc/ValProtonBunchTime.hh"
#include "Offline/Validation/inc/ValProtonBunchTimeMC.hh"
#include "Offline/Validation/inc/ValSimParticle.hh"
#include "Offline/Validation/inc/ValStatusG4.hh"
#include "Offline/Validation/inc/ValStepPointMC.hh"
#include "Offline/Validation/inc/ValSTMWaveformDigi.hh"
#include "Offline/Validation/inc/ValStrawDigi.hh"
#include "Offline/Validation/inc/ValStrawDigiADCWaveform.hh"
#include "Offline/Validation/inc/ValStrawDigiMC.hh"
#include "Offline/Validation/inc/ValStrawGasStep.hh"
#include "Offline/Validation/inc/ValStrawHit.hh"
#include "Offline/Validation/inc/ValStrawHitFlag.hh"
#include "Offline/Validation/inc/ValTimeCluster.hh"
#include "Offline/Validation/inc/ValTrackSummary.hh"
#include "Offline/Validation/inc/ValTriggerInfo.hh"
#include "Offline/Validation/inc/ValTriggerResults.hh"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/types/Atom.h"

namespace mu2e {

class Validation : public art::EDAnalyzer {
 public:
  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<int> validation_level{Name("validation_level"),
                                      Comment("validation level, 0 to 2"), 1};
  };

  // this line is required by art to allow the command line help print
  typedef art::EDAnalyzer::Table<Config> Parameters;

  explicit Validation(const Parameters& conf);
  void analyze(art::Event const& event) override;
  void beginJob() override;

 private:
  int _level;  // level=1 is a small number of histograms, 2 is more
  int _count;  // event count

  // ValXYZ are classes which contain a set of histograms for
  // validation of product XYZ.  They are in vectors, since we usually
  // have several instances of a product and we make histograms
  // for each instance.
  std::vector<std::shared_ptr<ValStatusG4>> _stat;
  std::vector<std::shared_ptr<ValProtonBunchIntensity>> _pbin;
  std::vector<std::shared_ptr<ValProtonBunchTime>> _pbtd;
  std::vector<std::shared_ptr<ValProtonBunchTimeMC>> _pbtm;
  std::vector<std::shared_ptr<ValEventWindowMarker>> _evwm;
  std::vector<std::shared_ptr<ValGenParticle>> _genp;
  std::vector<std::shared_ptr<ValSimParticle>> _simp;
  std::vector<std::shared_ptr<ValStepPointMC>> _spmc;
  std::vector<std::shared_ptr<ValCaloShowerStep>> _cals;
  std::vector<std::shared_ptr<ValCaloDigi>> _cald;
  std::vector<std::shared_ptr<ValCaloRecoDigi>> _calr;
  std::vector<std::shared_ptr<ValCaloHit>> _calh;
  std::vector<std::shared_ptr<ValCaloCluster>> _ccls;
  std::vector<std::shared_ptr<ValCrvStep>> _cvst;
  std::vector<std::shared_ptr<ValCrvDigi>> _cvdg;
  std::vector<std::shared_ptr<ValCrvDigiMC>> _cmdg;
  std::vector<std::shared_ptr<ValCrvRecoPulse>> _cvrp;
  std::vector<std::shared_ptr<ValCrvCoincidenceCluster>> _cvcc;
  std::vector<std::shared_ptr<ValSTMWaveformDigi>> _stmw;
  std::vector<std::shared_ptr<ValStrawGasStep>> _stgs;
  std::vector<std::shared_ptr<ValStrawDigi>> _stdg;
  std::vector<std::shared_ptr<ValStrawDigiADCWaveform>> _stdw;
  std::vector<std::shared_ptr<ValStrawDigiMC>> _stdm;
  std::vector<std::shared_ptr<ValStrawHit>> _stwh;
  std::vector<std::shared_ptr<ValBkgCluster>> _bgcl;
  std::vector<std::shared_ptr<ValBkgQual>> _bgql;
  std::vector<std::shared_ptr<ValTrackSummary>> _trks;
  std::vector<std::shared_ptr<ValHelixSeed>> _hxsd;
  std::vector<std::shared_ptr<ValKalSeed>> _klsd;
  std::vector<std::shared_ptr<ValStrawHitFlag>> _shfl;
  std::vector<std::shared_ptr<ValTimeCluster>> _tmcl;
  std::vector<std::shared_ptr<ValComboHit>> _stht;
  std::vector<std::shared_ptr<ValTriggerResults>> _trrs;
  std::vector<std::shared_ptr<ValTriggerInfo>> _tris;

  // Loop over the products of type T and
  // call fill() on validation histogram class V to make histograms.
  // It will create the root file subdirectories and
  // histograms on the first event
  template <class T, class V>
  int analyzeProduct(std::vector<std::shared_ptr<V>>& list,
                     art::Event const& event);

  art::ServiceHandle<art::TFileService> _tfs;
};

}  // namespace mu2e

mu2e::Validation::Validation(const Parameters& conf) :
    art::EDAnalyzer(conf), _level(conf().validation_level()), _count(0) {}

void mu2e::Validation::beginJob() {}

void mu2e::Validation::analyze(art::Event const& event) {
  analyzeProduct<StatusG4, ValStatusG4>(_stat, event);
  analyzeProduct<ProtonBunchIntensity, ValProtonBunchIntensity>(_pbin, event);
  analyzeProduct<ProtonBunchTime, ValProtonBunchTime>(_pbtd, event);
  analyzeProduct<ProtonBunchTimeMC, ValProtonBunchTimeMC>(_pbtm, event);
  analyzeProduct<EventWindowMarker, ValEventWindowMarker>(_evwm, event);
  analyzeProduct<GenParticleCollection, ValGenParticle>(_genp, event);
  analyzeProduct<SimParticleCollection, ValSimParticle>(_simp, event);
  analyzeProduct<StepPointMCCollection, ValStepPointMC>(_spmc, event);
  analyzeProduct<CaloShowerStepCollection, ValCaloShowerStep>(_cals, event);
  analyzeProduct<CaloDigiCollection, ValCaloDigi>(_cald, event);
  analyzeProduct<CaloRecoDigiCollection, ValCaloRecoDigi>(_calr, event);
  analyzeProduct<CaloHitCollection, ValCaloHit>(_calh, event);
  analyzeProduct<CaloClusterCollection, ValCaloCluster>(_ccls, event);
  analyzeProduct<CrvStepCollection, ValCrvStep>(_cvst, event);
  analyzeProduct<CrvDigiCollection, ValCrvDigi>(_cvdg, event);
  analyzeProduct<CrvDigiMCCollection, ValCrvDigiMC>(_cmdg, event);
  analyzeProduct<CrvRecoPulseCollection, ValCrvRecoPulse>(_cvrp, event);
  analyzeProduct<CrvCoincidenceClusterCollection, ValCrvCoincidenceCluster>(
      _cvcc, event);
  analyzeProduct<StrawGasStepCollection, ValStrawGasStep>(_stgs, event);
  analyzeProduct<StrawDigiCollection, ValStrawDigi>(_stdg, event);
  analyzeProduct<StrawDigiADCWaveformCollection, ValStrawDigiADCWaveform>(
      _stdw, event);
  analyzeProduct<StrawDigiMCCollection, ValStrawDigiMC>(_stdm, event);
  analyzeProduct<StrawHitCollection, ValStrawHit>(_stwh, event);
  analyzeProduct<StrawHitFlagCollection, ValStrawHitFlag>(_shfl, event);
  analyzeProduct<BkgClusterCollection, ValBkgCluster>(_bgcl, event);
  analyzeProduct<BkgQualCollection, ValBkgQual>(_bgql, event);
  analyzeProduct<ComboHitCollection, ValComboHit>(_stht, event);
  analyzeProduct<TimeClusterCollection, ValTimeCluster>(_tmcl, event);
  analyzeProduct<HelixSeedCollection, ValHelixSeed>(_hxsd, event);
  analyzeProduct<KalSeedCollection, ValKalSeed>(_klsd, event);
  analyzeProduct<TrackSummaryCollection, ValTrackSummary>(_trks, event);
  analyzeProduct<STMWaveformDigiCollection, ValSTMWaveformDigi>(_stmw, event);
  analyzeProduct<art::TriggerResults, ValTriggerResults>(_trrs, event);
  analyzeProduct<TriggerInfo, ValTriggerInfo>(_tris, event);
}

// make histograms for product type T by calling fill() on
// histogram class V
template <class T, class V>
int mu2e::Validation::analyzeProduct(std::vector<std::shared_ptr<V>>& list,
                                     art::Event const& event) {
  // get all instances of products of type T
  std::vector<art::Handle<T>> vah = event.getMany<T>();

  std::string name;
  // loop over the list of instances of products of this type
  for (auto const& ah : vah) {
    const art::Provenance* prov = ah.provenance();

    // the name of the root file directory holding these histos
    // is the className_moduleName_InstanceName for the instance
    std::string fcn = prov->friendlyClassName();
    // clean up some complicated "friendly" names
    if (fcn == "mu2e::SimParticleart::Ptrmu2e::MCTrajectorystd::map")
      fcn = "MCTrajectory";
    if (fcn == "mu2e::StrawHitFlagDetailmu2e::BitMaps") fcn = "StrawHitFlag";
    if (fcn == "mu2e::TrkQualDetailmu2e::MVAStructs") fcn = "TrkQual";
    if (fcn == "mu2e::BkgQualDetailmu2e::MVAStructs") fcn = "BkgQual";
    if (fcn.find("mu2e::", 0) == 0) fcn.erase(0, 6);
    if (fcn.find("art::", 0) == 0) fcn.erase(0, 5);

    std::string inst = prov->productInstanceName();
    if (inst.size() == 0) inst = "noName";
    // this verison has processname
    // name = fcn+"_"+prov->moduleLabel()+"_"+prov->processName()+"_"+ inst;
    name = fcn + "_" + prov->moduleLabel() + "_" + inst;

    // there is a TriggerResults for each procname, need to distinguish them
    if (fcn.find("TriggerResults", 0) == 0) {
      name = name + "_" + prov->processName();
    }

    // see if this instance of this product is already in our list
    // of products being histogrammed.  If not, add it to the list
    std::shared_ptr<V> prd = nullptr;
    for (auto ptr : list) {
      // if this instance of this product found in the list
      if (ptr->name().compare(name) == 0) prd = ptr;
    }
    // if not in the list, create a new set of histograms
    // for this product
    if (prd == nullptr) {
      prd = std::make_shared<V>(name);
      // create root file subdirectory
      const art::TFileDirectory& tfdir = _tfs->mkdir(name);
      // create histograms
      prd->declare(tfdir);
      // add it to the list of products being histogrammed
      list.push_back(prd);
    }
    // histogram this event for this product instance
    prd->fill(*ah, event);
  }
  return 0;
}

DEFINE_ART_MODULE(mu2e::Validation)
