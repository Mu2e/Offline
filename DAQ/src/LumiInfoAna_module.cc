// ======================================================================
//
// LumiInfoAna_module.cc : Analyze lumi stream info
//
// ======================================================================

// ROOT includes
#include "TH1.h"
#include "TH2.h"

// art includes
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "artdaq-core-mu2e/Data/EventHeader.hh"
#include "fhiclcpp/ParameterSet.h"

// Mu2e includes
#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"

// std includes
#include <iostream>
#include <memory>
#include <string>
#include <map>

// ======================================================================
namespace mu2e {

class LumiInfoAna : public art::EDAnalyzer {

public:
  struct Config {
    fhicl::Atom<int> diagLevel{fhicl::Name("diagLevel"), fhicl::Comment("diagnostic level"), 0};
  };

  explicit LumiInfoAna(const art::EDAnalyzer::Table<Config>& config);
  virtual ~LumiInfoAna() {}

  virtual void beginJob() override;
  virtual void endJob() override;
  virtual void endSubRun(const art::SubRun& sr) override;

  virtual void analyze(const art::Event& e) override;

  void fillCalo(std::vector<art::Handle<mu2e::IntensityInfosCalo>>& handles);
  void fillTrackerHits(std::vector<art::Handle<mu2e::IntensityInfosTrackerHits>>& handles);
  void fillTimeClusters(std::vector<art::Handle<mu2e::IntensityInfosTimeCluster>>& handles);
  void fillEventHeaders(std::vector<art::Handle<mu2e::EventHeaders>>& handles);

private:

  int _diagLevel;

  size_t _nevents;
  size_t _nsubrun;
  std::map<std::string, size_t> _counter_by_object;

  TH1* _hNCaloHits;
  TH1* _hCaloEnergy;
  TH1* _hCaphriEnergy;
  TH1* _hNCaphriHits;
  TH1* _hNTimeClusters;
  TH1* _hNTrackerHits;
  TH1* _hSubRuns;
};

// ======================================================================

LumiInfoAna::LumiInfoAna(const art::EDAnalyzer::Table<Config>& config) : art::EDAnalyzer{config}
    , _diagLevel(config().diagLevel())
    , _nevents(0)
    , _nsubrun(0)
{
}

//--------------------------------------------------------------------------------
// create the histograms
//--------------------------------------------------------------------------------
void LumiInfoAna::beginJob() {
  art::ServiceHandle<art::TFileService> tfs;
  art::TFileDirectory infoDir = tfs->mkdir("info");

  _hNCaloHits     = infoDir.make<TH1F>("hNCaloHits", "N(Calorimeter hits);N(hits);Events", 100, 0., 1000.);
  _hCaloEnergy    = infoDir.make<TH1F>("hCaloEnergy", "Calorimeter energy;Energy (MeV);Events", 100, 0., 6000.);
  _hCaphriEnergy  = infoDir.make<TH1F>("hCaphriEnergy", "CAPHRI energy;Energy (MeV);Events", 100, 0., 10.);
  _hNCaphriHits   = infoDir.make<TH1F>("hNCaphriHits", "N(CAPHRI hits);N(hits);Events", 100, 0., 100.);
  _hNTrackerHits  = infoDir.make<TH1F>("hNTrackerHits", "N(tracker hits);N(hits);Events", 100, 0., 10000.);
  _hNTimeClusters = infoDir.make<TH1F>("hNTimeClusters", "N(time clusters);N(time clusters);Events", 100, 0., 100.);
  _hSubRuns       = infoDir.make<TH1F>("hSubRuns", "Sub-Runs;Sub-Run;Events", 100, 0., 100.);
}

void LumiInfoAna::endJob() {
  std::cout << "LumiInfoAna::" << __func__ << ": Read information for an accumulated " << _nevents << " events\n";
}

//--------------------------------------------------------------------------------
void LumiInfoAna::fillCalo(std::vector<art::Handle<mu2e::IntensityInfosCalo>>& handles) {
  if(_diagLevel > 0) {
    std::cout << "LumiInfoAna::" << __func__ << ": Processing " << handles.size() << " handles"
              << std::endl;
  }
  for (const auto& handle : handles) {
    if (!handle.isValid() || handle->empty()) {
      if(_diagLevel > 0) {
        std::cout << "LumiInfoAna::" << __func__ << ": Invalid or empty calo info handle!"
                  << " Valid = " << handle.isValid() << " and size = " << ((handle.isValid()) ? ((int) handle->size()) : -1)
                  << std::endl;
      }
      continue;
    }
    _counter_by_object["calo"] += handle->size();
    for(auto& info : *handle) {
      _hNCaloHits  ->Fill(info.nCaloHits());
      _hCaloEnergy ->Fill(info.caloEnergy());
      _hNCaphriHits->Fill(info.nCaphriHits());
      for(size_t ihit = 0; ihit < info.nCaphriHits(); ++ihit) {
        double energy;
        int id;
        info.getCaphriHit(ihit, energy, id);
        _hCaphriEnergy->Fill(energy);
      }
    }
  }
}

//--------------------------------------------------------------------------------
void LumiInfoAna::fillTrackerHits(std::vector<art::Handle<mu2e::IntensityInfosTrackerHits>>& handles) {
  if(_diagLevel > 0) {
    std::cout << "LumiInfoAna::" << __func__ << ": Processing " << handles.size() << " handles"
              << std::endl;
  }
  for (const auto& handle : handles) {
    if (!handle.isValid() || handle->empty()) {
      if(_diagLevel > 0) {
        std::cout << "LumiInfoAna::" << __func__ << ": Invalid or empty tracker hit info handle!"
                  << " Valid = " << handle.isValid() << " and size = " << ((handle.isValid()) ? ((int) handle->size()) : -1)
                  << std::endl;
      }
      continue;
    }
    _counter_by_object["tracker"] += handle->size();
    for(auto& info : *handle) {
      _hNTrackerHits->Fill(info.nTrackerHits());
    }
  }
}

//--------------------------------------------------------------------------------
void LumiInfoAna::fillTimeClusters(std::vector<art::Handle<mu2e::IntensityInfosTimeCluster>>& handles) {
  if(_diagLevel > 0) {
    std::cout << "LumiInfoAna::" << __func__ << ": Processing " << handles.size() << " handles"
              << std::endl;
  }
  for (const auto& handle : handles) {
    if (!handle.isValid() || handle->empty()) {
      if(_diagLevel > 0) {
        std::cout << "LumiInfoAna::" << __func__ << ": Invalid or empty time clusters info handle!"
                  << " Valid = " << handle.isValid() << " and size = " << ((handle.isValid()) ? ((int) handle->size()) : -1)
                  << std::endl;
      }
      continue;
    }
    _counter_by_object["time_clusters"] += handle->size();
    for(auto& info : *handle) {
      _hNTimeClusters->Fill(info.nProtonTCs());
    }
  }
}

//--------------------------------------------------------------------------------
void LumiInfoAna::fillEventHeaders(std::vector<art::Handle<mu2e::EventHeaders>>& handles) {
  if(_diagLevel > 0) {
    std::cout << "LumiInfoAna::" << __func__ << ": Processing " << handles.size() << " handles"
              << std::endl;
  }
  for (const auto& handle : handles) {
    if (!handle.isValid() || handle->empty()) {
      if(_diagLevel > 0) {
        std::cout << "LumiInfoAna::" << __func__ << ": Invalid or empty event headers handle!"
                  << " Valid = " << handle.isValid() << " and size = " << ((handle.isValid()) ? ((int) handle->size()) : -1)
                  << std::endl;
      }
      continue;
    }
    // Count the number of events using the number of event headers found
    _nsubrun += handle->size();
    _counter_by_object["headers"] += handle->size();
    // Not yet histogramming any event header information
  }
}

//--------------------------------------------------------------------------------
void LumiInfoAna::endSubRun(const art::SubRun& sr) {
  art::SubRunNumber_t subrunNumber = sr.subRun();
  art::RunNumber_t    runNumber    = sr.run();

  //---------------------------------------
  // Retrieve the data

  // Calorimeter data
  auto caloHandles = sr.getMany<mu2e::IntensityInfosCalo>();
  // Tracker hit data
  auto trkHandles = sr.getMany<mu2e::IntensityInfosTrackerHits>();
  // Time cluster data
  auto tcHandles = sr.getMany<mu2e::IntensityInfosTimeCluster>();
  // Event header data
  auto headerHandles = sr.getMany<mu2e::EventHeaders>();
  if(_diagLevel > 0) {
    std::cout << "LumiInfoAna::" << __func__ << ": Subrun "
              << runNumber << ":" << subrunNumber
              << ": Retrieved " << caloHandles.size() << " calo handles"
              << " " << trkHandles.size() << " tracker hits handles"
              << " " << tcHandles.size() << " time cluster handles"
              << " " << headerHandles.size() << " header handles"
              << std::endl;
  }
  fillCalo(caloHandles);
  fillTrackerHits(trkHandles);
  fillTimeClusters(tcHandles);
  fillEventHeaders(headerHandles);
  _hSubRuns->Fill(Form("%i:%i", (int) runNumber, (int) subrunNumber), _nsubrun);
  _nevents += _nsubrun;
  if(_diagLevel > 0) {
    std::cout << "LumiInfoAna::" << __func__ << ": Subrun "
              << runNumber << ":" << subrunNumber
              << ": Subrun had " << _nsubrun << " events for a total of " << _nevents << " events"
              << std::endl;
    if(_diagLevel > 1) {
      std::cout << " Total counts by type:\n";
      for(auto entry : _counter_by_object) std::cout << "  " << entry.first.c_str()
                                                     << " : " << entry.second
                                                     << std::endl;
    }
  }
  _nsubrun = 0;
}

//--------------------------------------------------------------------------------
void LumiInfoAna::analyze(const art::Event& event) {
  // for printout use
  const art::EventNumber_t  eventNumber  = event.event();
  const art::SubRunNumber_t subrunNumber = event.subRun();
  const art::RunNumber_t    runNumber    = event.run();

  //---------------------------------------
  // Retrieve the data

  // Calorimeter data
  auto caloHandles = event.getMany<mu2e::IntensityInfosCalo>();
  // Tracker hit data
  auto trkHandles = event.getMany<mu2e::IntensityInfosTrackerHits>();
  // Time cluster data
  auto tcHandles = event.getMany<mu2e::IntensityInfosTimeCluster>();
  // Event header data
  auto headerHandles = event.getMany<mu2e::EventHeaders>();
  if(_diagLevel > 0) {
    std::cout << "LumiInfoAna::" << __func__ << ": Subrun "
              << runNumber << ":" << subrunNumber << ":" << eventNumber
              << ": Retrieved " << caloHandles.size() << " calo handles"
              << " " << trkHandles.size() << " tracker hits handles"
              << " " << tcHandles.size() << " time cluster handles"
              << " " << headerHandles.size() << " header handles"
              << std::endl;
  }
  fillCalo(caloHandles);
  fillTrackerHits(trkHandles);
  fillTimeClusters(tcHandles);
  fillEventHeaders(headerHandles);
}

} // end namespace mu2e

DEFINE_ART_MODULE(mu2e::LumiInfoAna)
