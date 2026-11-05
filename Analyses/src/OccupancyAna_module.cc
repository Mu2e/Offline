// OccupancyAna: simple per-event occupancy counting of tracker, calorimeter and CRV collections
// Writes one TTree entry per event with the size of each configured collection, plus per-plane
// (tracker) and per-disk (calorimeter) breakdowns. Needs no geometry, conditions or MC truth.
// Original author: Michael MacKenzie, 2026

// framework
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/OptionalAtom.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"

// Offline
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
#include "Offline/DataProducts/inc/CrystalId.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"

// ROOT
#include "TH1.h"
#include "TTree.h"

// c++
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace mu2e {

  class OccupancyAna : public art::EDAnalyzer {
  public:
    struct Config {
      using Name = fhicl::Name;
      using Comment = fhicl::Comment;
      // An unconfigured collection is skipped; its count is then -1 in the tree
      fhicl::OptionalAtom<art::InputTag> strawDigiTag   { Name("StrawDigiCollection")            , Comment("StrawDigi collection") };
      fhicl::OptionalAtom<art::InputTag> strawHitTag    { Name("StrawHitCollection")             , Comment("Single-straw ComboHit collection, e.g. makeSH") };
      fhicl::OptionalAtom<art::InputTag> comboHitTag    { Name("ComboHitCollection")             , Comment("Panel ComboHit collection, e.g. makePH") };
      fhicl::OptionalAtom<art::InputTag> caloDigiTag    { Name("CaloDigiCollection")             , Comment("CaloDigi collection") };
      fhicl::OptionalAtom<art::InputTag> caloRecoDigiTag{ Name("CaloRecoDigiCollection")         , Comment("CaloRecoDigi collection") };
      fhicl::OptionalAtom<art::InputTag> caloHitTag     { Name("CaloHitCollection")              , Comment("CaloHit collection") };
      fhicl::OptionalAtom<art::InputTag> caloClusterTag { Name("CaloClusterCollection")          , Comment("CaloCluster collection") };
      fhicl::OptionalAtom<art::InputTag> crvDigiTag     { Name("CrvDigiCollection")              , Comment("CrvDigi collection") };
      fhicl::OptionalAtom<art::InputTag> crvPulseTag    { Name("CrvRecoPulseCollection")         , Comment("CrvRecoPulse collection") };
      fhicl::OptionalAtom<art::InputTag> crvClusterTag  { Name("CrvCoincidenceClusterCollection"), Comment("CrvCoincidenceCluster collection") };
      fhicl::Atom<double> hitTimeMin { Name("HitTimeMin")      , Comment("Minimum hit time for the windowed (*_sel) counts (ns)"), -1.e9 };
      fhicl::Atom<double> hitTimeMax { Name("HitTimeMax")      , Comment("Maximum hit time for the windowed (*_sel) counts (ns)"),  1.e9 };
      fhicl::Atom<double> caloHitEMin{ Name("CaloHitEnergyMin"), Comment("Minimum CaloHit energy for the windowed count (MeV)")  ,  0.   };
      fhicl::Atom<bool>   requireAll { Name("RequireAll")      , Comment("Throw if a configured collection is missing (else count = -1)"), false };
      fhicl::Atom<int>    debugLevel { Name("DebugLevel")      , Comment("Debug level"), 0 };
    };

    // One counted collection
    struct Counter_t {
      std::string    name;
      art::InputTag  tag;
      bool           use = false; // configured
      int*           branch = nullptr; // tree branch the count is written to
      TH1*           hist = nullptr;
      long           nfound = 0; // events with the collection present
      long           nmissing = 0;
      double         sum = 0.;
      double         sum2 = 0.;
      int            max = 0;
    };

    // Output tree information
    // The per-plane and per-disk arrays are -1 when their collection is skipped or missing
    constexpr static int kNPlanes = StrawId::_nplanes;
    constexpr static int kNDisks  = CaloConst::_nDisk;
    constexpr static int kNQuickLookBins = 100; // initial binning, the axis extends as needed
    struct Tree_t {
      int run;
      int subrun;
      int event;
      int nstrawdigis;
      int nstrawhits;
      int ncombohits;
      int ncombohit_straws; // sum of the straw hits in the combo hits
      int ncalodigis;
      int ncalorecodigis;
      int ncalohits;
      int ncaloclusters;
      int ncrvdigis;
      int ncrvpulses;
      int ncrvclusters;
      int nstrawhits_sel; // in the time window
      int ncombohits_sel; // in the time window
      int ncalohits_sel;  // in the time window and above the energy threshold
      std::array<int, kNPlanes> strawdigis_plane;
      std::array<int, kNPlanes> strawhits_plane;
      std::array<int, kNDisks > calodigis_disk;
      std::array<int, kNDisks > calohits_disk;
    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit OccupancyAna(const Parameters& conf);
    void beginJob() override;
    void analyze(const art::Event& event) override;
    void endJob() override;

  private:
    // Retrieve a collection, recording its size; returns nullptr if skipped or missing
    template <class T> const T* getCollection(const art::Event& event, Counter_t& counter);
    static void setCounter(Counter_t& counter, const char* name, const fhicl::OptionalAtom<art::InputTag>& tag, int* branch);

    Counter_t strawDigis_;
    Counter_t strawHits_;
    Counter_t comboHits_;
    Counter_t caloDigis_;
    Counter_t caloRecoDigis_;
    Counter_t caloHits_;
    Counter_t caloClusters_;
    Counter_t crvDigis_;
    Counter_t crvPulses_;
    Counter_t crvClusters_;
    std::vector<Counter_t*> counters_;

    double hitTimeMin_;
    double hitTimeMax_;
    double caloHitEMin_;
    bool   requireAll_;
    int    debugLevel_;

    long   nevents_ = 0;
    Tree_t data_{};
    TTree* tree_ = nullptr;
  };

  //------------------------------------------------------------------------------------------------------------
  void OccupancyAna::setCounter(Counter_t& counter, const char* name, const fhicl::OptionalAtom<art::InputTag>& tag, int* branch) {
    counter.name   = name;
    counter.use    = tag(counter.tag);
    counter.branch = branch;
  }

  //------------------------------------------------------------------------------------------------------------
  OccupancyAna::OccupancyAna(const Parameters& conf)
    : art::EDAnalyzer(conf)
    , hitTimeMin_ (conf().hitTimeMin())
    , hitTimeMax_ (conf().hitTimeMax())
    , caloHitEMin_(conf().caloHitEMin())
    , requireAll_ (conf().requireAll())
    , debugLevel_ (conf().debugLevel())
  {
    setCounter(strawDigis_   , "strawdigis"   , conf().strawDigiTag   , &data_.nstrawdigis   );
    setCounter(strawHits_    , "strawhits"    , conf().strawHitTag    , &data_.nstrawhits    );
    setCounter(comboHits_    , "combohits"    , conf().comboHitTag    , &data_.ncombohits    );
    setCounter(caloDigis_    , "calodigis"    , conf().caloDigiTag    , &data_.ncalodigis    );
    setCounter(caloRecoDigis_, "calorecodigis", conf().caloRecoDigiTag, &data_.ncalorecodigis);
    setCounter(caloHits_     , "calohits"     , conf().caloHitTag     , &data_.ncalohits     );
    setCounter(caloClusters_ , "caloclusters" , conf().caloClusterTag , &data_.ncaloclusters );
    setCounter(crvDigis_     , "crvdigis"     , conf().crvDigiTag     , &data_.ncrvdigis     );
    setCounter(crvPulses_    , "crvpulses"    , conf().crvPulseTag    , &data_.ncrvpulses    );
    setCounter(crvClusters_  , "crvclusters"  , conf().crvClusterTag  , &data_.ncrvclusters  );
    counters_ = {&strawDigis_, &strawHits_, &comboHits_, &caloDigis_, &caloRecoDigis_, &caloHits_,
                 &caloClusters_, &crvDigis_, &crvPulses_, &crvClusters_};

    if(std::none_of(counters_.begin(), counters_.end(), [](const Counter_t* c) { return c->use; })) {
      throw cet::exception("BADCONFIG") << "OccupancyAna: no collections configured\n";
    }
  }

  //------------------------------------------------------------------------------------------------------------
  void OccupancyAna::beginJob() {
    art::ServiceHandle<art::TFileService> tfs;
    tree_ = tfs->make<TTree>("tree", "Per-event occupancy");
    tree_->Branch("run"             , &data_.run             );
    tree_->Branch("subrun"          , &data_.subrun          );
    tree_->Branch("event"           , &data_.event           );
    for(auto counter : counters_) tree_->Branch(("n" + counter->name).c_str(), counter->branch);
    tree_->Branch("ncombohit_straws", &data_.ncombohit_straws);
    tree_->Branch("nstrawhits_sel"  , &data_.nstrawhits_sel  );
    tree_->Branch("ncombohits_sel"  , &data_.ncombohits_sel  );
    tree_->Branch("ncalohits_sel"   , &data_.ncalohits_sel   );
    tree_->Branch("strawdigis_plane", data_.strawdigis_plane.data(), Form("strawdigis_plane[%i]/I", kNPlanes));
    tree_->Branch("strawhits_plane" , data_.strawhits_plane .data(), Form("strawhits_plane[%i]/I" , kNPlanes));
    tree_->Branch("calodigis_disk"  , data_.calodigis_disk  .data(), Form("calodigis_disk[%i]/I"  , kNDisks ));
    tree_->Branch("calohits_disk"   , data_.calohits_disk   .data(), Form("calohits_disk[%i]/I"   , kNDisks ));

    // Quick-look collection size histograms, for the configured collections only
    for(auto counter : counters_) {
      if(!counter->use) continue;
      counter->hist = tfs->make<TH1D>(("n" + counter->name).c_str(),
                                      Form("N(%s) per event;N(%s);Events", counter->tag.encode().c_str(), counter->name.c_str()),
                                      kNQuickLookBins, 0., kNQuickLookBins);
      counter->hist->SetCanExtend(TH1::kXaxis);
    }
  }

  //------------------------------------------------------------------------------------------------------------
  template <class T> const T* OccupancyAna::getCollection(const art::Event& event, Counter_t& counter) {
    *counter.branch = -1;
    if(!counter.use) return nullptr;
    art::Handle<T> handle;
    event.getByLabel(counter.tag, handle);
    if(!handle.isValid()) {
      if(requireAll_) throw cet::exception("BADINPUT") << "OccupancyAna: missing " << counter.name
                                                       << " collection with tag " << counter.tag.encode() << "\n";
      if(counter.nmissing == 0) printf("[OccupancyAna::%s] Warning: %s collection %s not found, its count is set to -1\n",
                                       __func__, counter.name.c_str(), counter.tag.encode().c_str());
      ++counter.nmissing;
      return nullptr;
    }
    const int n = handle->size();
    *counter.branch = n;
    ++counter.nfound;
    counter.sum  += n;
    counter.sum2 += double(n)*n;
    counter.max   = std::max(counter.max, n);
    counter.hist->Fill(n);
    return handle.product();
  }

  //------------------------------------------------------------------------------------------------------------
  void OccupancyAna::analyze(const art::Event& event) {
    ++nevents_;
    data_.run    = static_cast<int>(event.run());
    data_.subrun = static_cast<int>(event.subRun());
    data_.event  = static_cast<int>(event.event());
    data_.ncombohit_straws = -1;
    data_.nstrawhits_sel   = -1;
    data_.ncombohits_sel   = -1;
    data_.ncalohits_sel    = -1;
    data_.strawdigis_plane.fill(-1);
    data_.strawhits_plane .fill(-1);
    data_.calodigis_disk  .fill(-1);
    data_.calohits_disk   .fill(-1);

    const auto inWindow = [&](const double t) { return t >= hitTimeMin_ && t <= hitTimeMax_; };

    // Tracker
    if(const auto col = getCollection<StrawDigiCollection>(event, strawDigis_)) {
      data_.strawdigis_plane.fill(0);
      for(const auto& digi : *col) {
        const int plane = digi.strawId().plane();
        if(plane >= 0 && plane < kNPlanes) ++data_.strawdigis_plane[plane];
      }
    }
    if(const auto col = getCollection<ComboHitCollection>(event, strawHits_)) {
      data_.strawhits_plane.fill(0);
      data_.nstrawhits_sel = 0;
      for(const auto& hit : *col) {
        const int plane = hit.strawId().plane();
        if(plane >= 0 && plane < kNPlanes) ++data_.strawhits_plane[plane];
        if(inWindow(hit.correctedTime())) ++data_.nstrawhits_sel;
      }
    }
    if(const auto col = getCollection<ComboHitCollection>(event, comboHits_)) {
      data_.ncombohit_straws = 0;
      data_.ncombohits_sel   = 0;
      for(const auto& hit : *col) {
        data_.ncombohit_straws += hit.nStrawHits();
        if(inWindow(hit.correctedTime())) ++data_.ncombohits_sel;
      }
    }

    // Calorimeter
    if(const auto col = getCollection<CaloDigiCollection>(event, caloDigis_)) {
      data_.calodigis_disk.fill(0);
      for(const auto& digi : *col) {
        const CaloSiPMId id(digi.SiPMID());
        if(id.isCrystal()) ++data_.calodigis_disk[id.disk()];
      }
    }
    getCollection<CaloRecoDigiCollection>(event, caloRecoDigis_);
    if(const auto col = getCollection<CaloHitCollection>(event, caloHits_)) {
      data_.calohits_disk.fill(0);
      data_.ncalohits_sel = 0;
      for(const auto& hit : *col) {
        const CrystalId id(hit.crystalID());
        if(id.isValid()) ++data_.calohits_disk[id.disk()];
        if(inWindow(hit.time()) && hit.energyDep() >= caloHitEMin_) ++data_.ncalohits_sel;
      }
    }
    getCollection<CaloClusterCollection>(event, caloClusters_);

    // CRV
    getCollection<CrvDigiCollection>(event, crvDigis_);
    getCollection<CrvRecoPulseCollection>(event, crvPulses_);
    getCollection<CrvCoincidenceClusterCollection>(event, crvClusters_);

    if(debugLevel_ > 0) {
      printf("[OccupancyAna::%s] Run %i subrun %i event %i:", __func__, data_.run, data_.subrun, data_.event);
      for(auto counter : counters_) if(counter->use) printf(" N(%s) = %i", counter->name.c_str(), *counter->branch);
      printf("\n");
    }

    tree_->Fill();
  }

  //------------------------------------------------------------------------------------------------------------
  void OccupancyAna::endJob() {
    // Summary table, one line per configured collection
    printf("[OccupancyAna::%s] Occupancy summary over %li events:\n", __func__, nevents_);
    printf("  %-14s %-28s %8s %8s %12s %12s %8s\n", "collection", "tag", "found", "missing", "mean", "rms", "max");
    for(auto counter : counters_) {
      if(!counter->use) continue;
      const double nfound = static_cast<double>(counter->nfound);
      const double mean = (counter->nfound > 0) ? counter->sum/nfound : 0.;
      const double rms  = (counter->nfound > 0) ? std::sqrt(std::max(0., counter->sum2/nfound - mean*mean)) : 0.;
      printf("  %-14s %-28s %8li %8li %12.3f %12.3f %8i\n", counter->name.c_str(), counter->tag.encode().c_str(),
             counter->nfound, counter->nmissing, mean, rms, counter->max);
    }
  }
}

DEFINE_ART_MODULE(mu2e::OccupancyAna)
