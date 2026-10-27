// OccupancyAna: simple per-event occupancy counting of tracker, calorimeter and CRV collections
// Writes one TTree entry per event with the size of each configured collection, plus per-plane
// (tracker) and per-disk (calorimeter) breakdowns. Needs no geometry, conditions or MC truth.
// Original author: Michael MacKenzie, 2026

// framework
#include "fhiclcpp/types/Atom.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art_root_io/TFileService.h"
#include "cetlib_except/exception.h"

// Offline
#include "Offline/DataProducts/inc/CaloConst.hh"
#include "Offline/DataProducts/inc/CaloSiPMId.hh"
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
      // An empty tag skips that collection; its count is then -1 in the tree
      fhicl::Atom<std::string> strawDigiTag   { Name("StrawDigiCollection")            , Comment("StrawDigi collection (empty to skip)")            , "" };
      fhicl::Atom<std::string> strawHitTag    { Name("StrawHitCollection")             , Comment("Single-straw ComboHit collection, e.g. makeSH (empty to skip)"), "" };
      fhicl::Atom<std::string> comboHitTag    { Name("ComboHitCollection")             , Comment("Panel ComboHit collection, e.g. makePH (empty to skip)"), "" };
      fhicl::Atom<std::string> caloDigiTag    { Name("CaloDigiCollection")             , Comment("CaloDigi collection (empty to skip)")             , "" };
      fhicl::Atom<std::string> caloRecoDigiTag{ Name("CaloRecoDigiCollection")         , Comment("CaloRecoDigi collection (empty to skip)")         , "" };
      fhicl::Atom<std::string> caloHitTag     { Name("CaloHitCollection")              , Comment("CaloHit collection (empty to skip)")              , "" };
      fhicl::Atom<std::string> caloClusterTag { Name("CaloClusterCollection")          , Comment("CaloCluster collection (empty to skip)")          , "" };
      fhicl::Atom<std::string> crvDigiTag     { Name("CrvDigiCollection")              , Comment("CrvDigi collection (empty to skip)")              , "" };
      fhicl::Atom<std::string> crvPulseTag    { Name("CrvRecoPulseCollection")         , Comment("CrvRecoPulse collection (empty to skip)")         , "" };
      fhicl::Atom<std::string> crvClusterTag  { Name("CrvCoincidenceClusterCollection"), Comment("CrvCoincidenceCluster collection (empty to skip)"), "" };
      fhicl::Atom<double>      hitTimeMin     { Name("HitTimeMin")                     , Comment("Minimum hit time for the windowed (*_sel) counts (ns)"), -1.e9 };
      fhicl::Atom<double>      hitTimeMax     { Name("HitTimeMax")                     , Comment("Maximum hit time for the windowed (*_sel) counts (ns)"),  1.e9 };
      fhicl::Atom<double>      caloHitEMin    { Name("CaloHitEnergyMin")               , Comment("Minimum CaloHit energy for the windowed count (MeV)")  ,  0.   };
      fhicl::Atom<bool>        requireAll     { Name("RequireAll")                     , Comment("Throw if a configured collection is missing (else count = -1)"), false };
      fhicl::Atom<int>         debugLevel     { Name("DebugLevel")                     , Comment("Debug level"), 0 };
    };

    // One counted collection
    struct Counter_t {
      std::string name;
      std::string tag;
      int*  branch = nullptr; // tree branch the count is written to
      TH1*  hist = nullptr;
      long  nfound = 0; // events with the collection present
      long  nmissing = 0;
      double sum = 0.;
      double sum2 = 0.;
      int   max = 0;
      bool  use() const { return !tag.empty(); }
    };

    // Output tree information
    enum {kNPlanes = StrawId::_nplanes, kNDisks = CaloConst::_nDisk};
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
      int strawdigis_plane[kNPlanes];
      int strawhits_plane [kNPlanes];
      int calodigis_disk  [kNDisks];
      int calohits_disk   [kNDisks];
    };

    using Parameters = art::EDAnalyzer::Table<Config>;
    explicit OccupancyAna(const Parameters& conf);
    void beginJob() override;
    void analyze(const art::Event& event) override;
    void endJob() override;

  private:
    // Retrieve a collection, recording its size; returns nullptr if skipped or missing
    template <class T> const T* getCollection(const art::Event& event, Counter_t& counter);

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
    Tree_t data_;
    TTree* tree_ = nullptr;
  };

  //------------------------------------------------------------------------------------------------------------
  OccupancyAna::OccupancyAna(const Parameters& conf)
    : art::EDAnalyzer(conf)
    , hitTimeMin_ (conf().hitTimeMin())
    , hitTimeMax_ (conf().hitTimeMax())
    , caloHitEMin_(conf().caloHitEMin())
    , requireAll_ (conf().requireAll())
    , debugLevel_ (conf().debugLevel())
  {
    strawDigis_   .name = "strawdigis"    ; strawDigis_   .tag = conf().strawDigiTag()   ; strawDigis_   .branch = &data_.nstrawdigis   ;
    strawHits_    .name = "strawhits"     ; strawHits_    .tag = conf().strawHitTag()    ; strawHits_    .branch = &data_.nstrawhits    ;
    comboHits_    .name = "combohits"     ; comboHits_    .tag = conf().comboHitTag()    ; comboHits_    .branch = &data_.ncombohits    ;
    caloDigis_    .name = "calodigis"     ; caloDigis_    .tag = conf().caloDigiTag()    ; caloDigis_    .branch = &data_.ncalodigis    ;
    caloRecoDigis_.name = "calorecodigis" ; caloRecoDigis_.tag = conf().caloRecoDigiTag(); caloRecoDigis_.branch = &data_.ncalorecodigis;
    caloHits_     .name = "calohits"      ; caloHits_     .tag = conf().caloHitTag()     ; caloHits_     .branch = &data_.ncalohits     ;
    caloClusters_ .name = "caloclusters"  ; caloClusters_ .tag = conf().caloClusterTag() ; caloClusters_ .branch = &data_.ncaloclusters ;
    crvDigis_     .name = "crvdigis"      ; crvDigis_     .tag = conf().crvDigiTag()     ; crvDigis_     .branch = &data_.ncrvdigis     ;
    crvPulses_    .name = "crvpulses"     ; crvPulses_    .tag = conf().crvPulseTag()    ; crvPulses_    .branch = &data_.ncrvpulses    ;
    crvClusters_  .name = "crvclusters"   ; crvClusters_  .tag = conf().crvClusterTag()  ; crvClusters_  .branch = &data_.ncrvclusters  ;
    counters_ = {&strawDigis_, &strawHits_, &comboHits_, &caloDigis_, &caloRecoDigis_, &caloHits_,
                 &caloClusters_, &crvDigis_, &crvPulses_, &crvClusters_};

    if(std::none_of(counters_.begin(), counters_.end(), [](const Counter_t* c) { return c->use(); })) {
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
    tree_->Branch("strawdigis_plane", data_.strawdigis_plane, Form("strawdigis_plane[%i]/I", int(kNPlanes)));
    tree_->Branch("strawhits_plane" , data_.strawhits_plane , Form("strawhits_plane[%i]/I" , int(kNPlanes)));
    tree_->Branch("calodigis_disk"  , data_.calodigis_disk  , Form("calodigis_disk[%i]/I"  , int(kNDisks )));
    tree_->Branch("calohits_disk"   , data_.calohits_disk   , Form("calohits_disk[%i]/I"   , int(kNDisks )));

    // Quick-look collection size histograms, for the configured collections only
    for(auto counter : counters_) {
      if(!counter->use()) continue;
      counter->hist = tfs->make<TH1D>(("n" + counter->name).c_str(),
                                      Form("N(%s) per event;N(%s);Events", counter->tag.c_str(), counter->name.c_str()),
                                      100, 0., 100.);
      counter->hist->SetCanExtend(TH1::kXaxis);
    }
  }

  //------------------------------------------------------------------------------------------------------------
  template <class T> const T* OccupancyAna::getCollection(const art::Event& event, Counter_t& counter) {
    *counter.branch = -1;
    if(!counter.use()) return nullptr;
    art::Handle<T> handle;
    event.getByLabel(art::InputTag(counter.tag), handle);
    if(!handle.isValid()) {
      if(requireAll_) throw cet::exception("BADINPUT") << "OccupancyAna: missing " << counter.name
                                                       << " collection with tag " << counter.tag << "\n";
      if(counter.nmissing == 0) printf("[OccupancyAna::%s] Warning: %s collection %s not found, its count is set to -1\n",
                                       __func__, counter.name.c_str(), counter.tag.c_str());
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
    data_.run    = event.run();
    data_.subrun = event.subRun();
    data_.event  = event.event();
    data_.ncombohit_straws = -1;
    data_.nstrawhits_sel   = -1;
    data_.ncombohits_sel   = -1;
    data_.ncalohits_sel    = -1;
    for(int i = 0; i < kNPlanes; ++i) { data_.strawdigis_plane[i] = 0; data_.strawhits_plane[i] = 0; }
    for(int i = 0; i < kNDisks ; ++i) { data_.calodigis_disk  [i] = 0; data_.calohits_disk  [i] = 0; }

    const auto inWindow = [&](const double t) { return t >= hitTimeMin_ && t <= hitTimeMax_; };

    // Tracker
    if(const auto col = getCollection<StrawDigiCollection>(event, strawDigis_)) {
      for(const auto& digi : *col) {
        const int plane = digi.strawId().plane();
        if(plane >= 0 && plane < kNPlanes) ++data_.strawdigis_plane[plane];
      }
    }
    if(const auto col = getCollection<ComboHitCollection>(event, strawHits_)) {
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
      for(const auto& digi : *col) {
        const CaloSiPMId id(digi.SiPMID());
        if(id.isCrystal()) ++data_.calodigis_disk[id.disk()];
      }
    }
    getCollection<CaloRecoDigiCollection>(event, caloRecoDigis_);
    if(const auto col = getCollection<CaloHitCollection>(event, caloHits_)) {
      data_.ncalohits_sel = 0;
      for(const auto& hit : *col) {
        const int disk = hit.crystalID() / CaloConst::_nCrystalPerDisk;
        if(disk >= 0 && disk < kNDisks) ++data_.calohits_disk[disk];
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
      for(auto counter : counters_) if(counter->use()) printf(" N(%s) = %i", counter->name.c_str(), *counter->branch);
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
      if(!counter->use()) continue;
      const double mean = (counter->nfound > 0) ? counter->sum/counter->nfound : 0.;
      const double rms  = (counter->nfound > 0) ? std::sqrt(std::max(0., counter->sum2/counter->nfound - mean*mean)) : 0.;
      printf("  %-14s %-28s %8li %8li %12.3f %12.3f %8i\n", counter->name.c_str(), counter->tag.c_str(),
             counter->nfound, counter->nmissing, mean, rms, counter->max);
    }
  }
}

DEFINE_ART_MODULE(mu2e::OccupancyAna)
