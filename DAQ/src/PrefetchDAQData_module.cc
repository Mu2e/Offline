// framework
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "fhiclcpp/ParameterSet.h"

// data
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"

// DAQ data
// #include "artdaq-core-mu2e/Overlays/FragmentType.hh"
// #include "artdaq-core-mu2e/Overlays/TrackerDataDecoder.hh"
#include <artdaq-core/Data/Fragment.hh>

namespace mu2e {

class PrefetchDAQData : public art::EDProducer {

protected:
public:
  explicit PrefetchDAQData(fhicl::ParameterSet const&);
  virtual ~PrefetchDAQData();
  virtual void beginJob();
  virtual void produce(art::Event& e);

  void fake_access(const CaloDigi& Digi);
  void fake_access(const StrawDigi& Digi);
  void fake_access(const artdaq::Fragment& Frag);

private:
  bool findData(const art::Event& e);
  // control flags
  int _debugLevel;
  int _fetchCaloDigis;
  int _fetchStrawDigis;
  int _fetchCaloFragments;
  int _fetchTrkFragments;
  int _fetchAllFragments;

  art::InputTag _cdTag;
  art::InputTag _sdTag;
  art::InputTag _cfTag;
  art::InputTag _tfTag;
  // cache of event objects
  const CaloDigiCollection* _cdcol;
  const StrawDigiCollection* _sdcol;
  const artdaq::Fragments* _cfcol;
  const artdaq::Fragments* _tfcol;
  int _eventNum;
};

//-----------------------------------------------------------------------------
PrefetchDAQData::PrefetchDAQData(fhicl::ParameterSet const& pset) :
    art::EDProducer(pset), _debugLevel(pset.get<int>("debugLevel")),
    _fetchCaloDigis(pset.get<int>("fetchCaloDigis")),
    _fetchStrawDigis(pset.get<int>("fetchStrawDigis")),
    _fetchCaloFragments(pset.get<int>("fetchCaloFragments")),
    _fetchTrkFragments(pset.get<int>("fetchTrkFragments")),
    _fetchAllFragments(pset.get<int>("fetchAllFragments")),
    _cdTag(pset.get<std::string>("caloDigiCollectionTag")),
    _sdTag(pset.get<std::string>("strawDigiCollectionTag")),
    _cfTag(pset.get<std::string>("caloFragmentTag")),
    _tfTag(pset.get<std::string>("trkFragmentTag")) {}

PrefetchDAQData::~PrefetchDAQData() {}

//-----------------------------------------------------------------------------
void PrefetchDAQData::beginJob() {}

//-----------------------------------------------------------------------------
void PrefetchDAQData::fake_access(const CaloDigi& Hit) {}

//-----------------------------------------------------------------------------
void PrefetchDAQData::fake_access(const StrawDigi& X) {}
//-----------------------------------------------------------------------------
void PrefetchDAQData::fake_access(const artdaq::Fragment& Frag) {}

//-----------------------------------------------------------------------------
bool PrefetchDAQData::findData(const art::Event& evt) {
  _cdcol = 0;
  _sdcol = 0;
  _cfcol = 0;
  _tfcol = 0;

  if (_fetchCaloDigis) {
    auto cdH = evt.getValidHandle<CaloDigiCollection>(_cdTag);
    _cdcol = cdH.product();
  }

  if (_fetchTrkFragments) {
    auto tfH = evt.getValidHandle<artdaq::Fragments>(_tfTag);
    _tfcol = tfH.product();
  }

  if (_fetchCaloFragments) {
    auto cfH = evt.getValidHandle<artdaq::Fragments>(_cfTag);
    _cfcol = cfH.product();
  }

  if (_fetchStrawDigis) {
    auto sdH = evt.getValidHandle<StrawDigiCollection>(_sdTag);
    _sdcol = sdH.product();
  }

  if (_fetchAllFragments) {
    auto allH = evt.getMany<std::vector<artdaq::Fragment>>();
    for (auto& fs : allH) {
      for (auto& frag : *fs) {
        fake_access(frag);
      }
    }
  }
  //-----------------------------------------------------------------------------
  // prefetch data
  //-----------------------------------------------------------------------------
  if (_cdcol) {
    int ncd = _cdcol->size();
    for (int i = 0; i < ncd; ++i) {
      const CaloDigi& cd = _cdcol->at(i);
      fake_access(cd);
    }
  }

  if (_sdcol) {
    int nsd = _sdcol->size();
    for (int i = 0; i < nsd; ++i) {
      const StrawDigi& sdigi = _sdcol->at(i);
      fake_access(sdigi);
    }
  }

  if (_cfcol) {
    int ncf = _cfcol->size();
    for (int i = 0; i < ncf; ++i) {
      auto frag = _cfcol->at(i);
      fake_access(frag);
    }
  }
  if (_tfcol) {
    int ntf = _tfcol->size();
    for (int i = 0; i < ntf; ++i) {
      auto frag = _tfcol->at(i);
      fake_access(frag);
    }
  }

  return true;
}

//-----------------------------------------------------------------------------
void PrefetchDAQData::produce(art::Event& Event) {

  _eventNum = Event.event();

  if (_debugLevel > 0)
    printf("[PrefetchDAQData::produce] event number: %10i\n", _eventNum);

  findData(Event);
}

//------------------------------------------------------------------------------
// Part of the magic that makes this class a module.
DEFINE_ART_MODULE(PrefetchDAQData)

} // namespace mu2e
