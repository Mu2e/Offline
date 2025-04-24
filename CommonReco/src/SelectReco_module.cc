//
//  Create a subset of digi and hit objects associated with given
//  reconstruction output objects. This also creates reco summary objects
//
// Original author: Dave Brown (LBNL) Feb 2025
//
// art
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "canvas/Utilities/InputTag.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Selector.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
// mu2e data products
#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/DataProducts/inc/VirtualDetectorId.hh"
#include "Offline/DataProducts/inc/IndexMap.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/RecoCount.hh"
#include "Offline/RecoDataProducts/inc/KKLoopHelix.hh"
#include "Offline/RecoDataProducts/inc/KKCentralHelix.hh"
#include "Offline/RecoDataProducts/inc/KKLine.hh"
// Utilities
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"
#include "Offline/TrackerGeom/inc/Straw.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/DetectorSystem.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
// C++
#include <vector>
#include <memory>
#include <iostream>
#include <set>
#include <string>
#include <algorithm>

namespace mu2e {
  class SelectReco : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
        fhicl::Atom<int>  debug                   { Name("debugLevel"),                     Comment("Debug Level"), 0};
        fhicl::Atom<art::InputTag> CCC            { Name("CaloClusterCollection"),          Comment("CaloClusterCollection")};
        fhicl::Atom<art::InputTag> CrvCCC         { Name("CrvCoincidenceClusterCollection"),Comment("CrvCoincidenceClusterCollections")};
        fhicl::Atom<art::InputTag> SDC            { Name("StrawDigiCollection"),            Comment("StrawDigiCollection")};
        fhicl::Atom<art::InputTag> CHC            { Name("ComboHitCollection"),             Comment("ComboHitCollection for the original StrawHits (not Panel hits)")};
        fhicl::Atom<art::InputTag> CDC            { Name("CaloDigiCollection"),             Comment("CaloDigiCollection")};
        fhicl::Atom<art::InputTag> CRVDC          { Name("CrvDigiCollection"),              Comment("CrvDigiCollection")};
        fhicl::Sequence<art::InputTag> KalSeeds   { Name("KalSeedCollections"),             Comment("KalSeedCollections")};
// selection parameters
        fhicl::Atom<double> CCME                  { Name("CaloClusterMinE"),                Comment("Minimum energy CaloCluster to save (MeV)")};
        fhicl::Atom<bool> saveNearby              { Name("SaveNearbyDigis"),                Comment("Save StrawDigis near the fit")};

      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit SelectReco(const Parameters& conf);
      void produce(art::Event& evt) override;

    private:
      typedef std::set<StrawDigiIndex> SDIS;
      // utility functions
      void fillStrawHitCounts(ComboHitCollection const& chc, RecoCount& nrec);
      void fillTrk           (art::Event& event, std::set<art::Ptr<CaloCluster> >& ccptrs, RecoCount& nrec);
      void fillCrv           (art::Event& event, RecoCount& nrec);
      void fillCalo          (art::Event& event, std::set<art::Ptr<CaloCluster> >& ccptrs, RecoCount& nrec);

      int debug_;
      art::InputTag ccct_, sdct_, chct_, cdct_, crvcct_, crvdct_;
      std::vector<art::InputTag> kscts_;
      double ccme_;
      bool saveunused_;
  };

  SelectReco::SelectReco(const Parameters& config )  :
    art::EDProducer{config},
    debug_(config().debug()),
    ccct_(config().CCC()),
    sdct_(config().SDC()),
    chct_(config().CHC()),
    cdct_(config().CDC()),
    crvcct_(config().CrvCCC()),
    crvdct_(config().CRVDC()),
    kscts_(config().KalSeeds()),
    ccme_(config().CCME()),
    saveunused_(config().saveNearby())
   {
      consumes<StrawDigiCollection>(sdct_);
      consumes<StrawDigiADCWaveformCollection>(sdct_);
      consumes<ComboHitCollection>(chct_);
      consumes<CaloDigiCollection>(cdct_);
      consumes<CaloClusterCollection>(ccct_);
      consumes<CrvDigiCollection>(crvdct_);
      consumesMany<KalSeedCollection>();
      consumes<CrvCoincidenceClusterCollection>(crvcct_);
      produces <IndexMap>("StrawDigiMap");
      produces <IndexMap>("CrvDigiMap");
      produces <CaloDigiCollection>();
      produces <CrvDigiCollection>();
      produces <CrvRecoPulseCollection>();
      produces <CrvCoincidenceClusterCollection>();
      produces <StrawDigiCollection>();
      produces <StrawDigiADCWaveformCollection>();
      produces <RecoCount>();

      if (debug_ > 0) {
        std::cout << "Using KalSeed collections from ";
        for (auto const& kff : kscts_) std::cout << kff << " " << std::endl;
      }
    }



  void SelectReco::produce(art::Event& event) {
    std::unique_ptr<RecoCount> nrec(new RecoCount);
    std::set<art::Ptr<CaloCluster> > ccptrs;
    fillTrk(event,ccptrs,*nrec.get());
    fillCrv(event, *nrec.get());
    fillCalo(event, ccptrs, *nrec.get());
    event.put(std::move(nrec));
  }

  void SelectReco::fillTrk( art::Event& event, std::set<art::Ptr<CaloCluster> >& ccptrs, RecoCount& nrec) {
   // Tracker-reated data products
    auto sdch = event.getValidHandle<StrawDigiCollection>(sdct_);
    auto const& sdc = *sdch;
    auto sdadcch = event.getValidHandle<StrawDigiADCWaveformCollection>(sdct_);
    auto const& sdadcc = *sdadcch;
    auto chch = event.getValidHandle<ComboHitCollection>(chct_);
    auto const& chc = *chch;
    // create products related to the reconstruction output
    std::unique_ptr<StrawDigiCollection> ssdc(new StrawDigiCollection);
    std::unique_ptr<StrawDigiADCWaveformCollection> ssdadcc(new StrawDigiADCWaveformCollection);
    // index maps between original collections and pruned collections
    std::unique_ptr<IndexMap> sdim(new IndexMap);
    // straw digi indices that are referenced by the tracks, or are 'close to' the track
    SDIS sdindices;
    // loop over input KalSeeds
    for (auto const& ksct : kscts_){
      auto ksch = event.getHandle<KalSeedCollection>(ksct);
      if(ksch.isValid()){
        auto const& ksc = *ksch;
        if(debug_ > 1) std::cout << "Found " << ksc.size() << " seeds from collection " << ksct << std::endl;
        for(auto const& ks : ksc) {
          // add Digis associated with the hits used in this fit
          for(auto const& hit : ks.hits() ) sdindices.insert(hit.index());
          // record the CaloCluster associated with this ks (if any)
          if(ks.hasCaloCluster())ccptrs.insert(ks.caloCluster());
        }
      }
    }
    // fill the StrawIndex map with the complete list of indices.
    StrawDigiIndex sdcount(0);
    ssdc->reserve(sdindices.size());
    ssdadcc->reserve(sdindices.size());
    for(auto sdindex : sdindices){
      sdim->addElement(sdindex,sdcount++);
      // deep-copy the selected StrawDigis
      ssdc->push_back(sdc[sdindex]);
      ssdadcc->push_back(sdadcc[sdindex]);
    }
    if(debug_ > 1) std::cout << "Selected " << sdcount << " StrawDigis" << std::endl;

    // fill detailed StrawHit counts
    fillStrawHitCounts(chc,nrec);
    event.put(std::move(sdim),"StrawDigiMap");
    event.put(std::move(ssdc));
    event.put(std::move(ssdadcc));
  }

  void SelectReco::fillStrawHitCounts(ComboHitCollection const& chc, RecoCount& nrec) {

    for(const auto& ch : chc) {
      auto const& shf = ch.flag();
      if(shf.hasAllProperties(StrawHitFlag::energysel))++nrec._nshfesel;
      if(shf.hasAllProperties(StrawHitFlag::radsel))++nrec._nshfrsel;
      if(shf.hasAllProperties(StrawHitFlag::timesel))++nrec._nshftsel;
      if(shf.hasAllProperties(StrawHitFlag::bkg))++nrec._nshfbkg;
      if(shf.hasAllProperties(StrawHitFlag::trksel))++nrec._nshftpk;
    }
    nrec._nstrawdigi = chc.size();
    // fill straw hit time histogram
    for(auto const& ch : chc)nrec._shthist.fill(ch.time());
  }
  void SelectReco::fillCrv(art::Event& event, RecoCount& nrec) {
    // find Crv data in event
    auto crvdch = event.getValidHandle<CrvDigiCollection>(crvdct_);
    auto const& crvdc = *crvdch;
    auto CrvRecoPulseCollectionPID = event.getProductID<CrvRecoPulseCollection>();
    auto CrvRecoPulseCollectionGetter = event.productGetter(CrvRecoPulseCollectionPID);
    // create new Crv collections
    std::unique_ptr<CrvDigiCollection> scrvdc(new CrvDigiCollection);
    std::unique_ptr<CrvRecoPulseCollection> scrvrpc(new CrvRecoPulseCollection);
    std::unique_ptr<CrvCoincidenceClusterCollection> scrvccc(new CrvCoincidenceClusterCollection);
    std::unique_ptr<IndexMap> crvdim(new IndexMap);

    std::set<uint16_t> crvindices;
    auto crvccch = event.getValidHandle<CrvCoincidenceClusterCollection>(crvcct_);
    auto const& crvccc = *crvccch;
    // loop over CrvCoincidenceClusters
    for(auto const& crvcc: crvccc) {
      std::vector<art::Ptr<CrvRecoPulse>> pulses;
      for(auto const& crvrp : crvcc.GetCrvRecoPulses()){
        // deep-copy the pulses used in coincidences: the digi indices are updated later
        // the map must be used to connect them
        scrvrpc->push_back(*crvrp);
        auto crvrpp = art::Ptr<CrvRecoPulse>(CrvRecoPulseCollectionPID,scrvrpc->size()-1,CrvRecoPulseCollectionGetter);
        pulses.push_back(crvrpp);
        for(auto index : crvrp->GetWaveformIndices()){
          crvindices.insert(index);
        }
      }
      // deep-copy the coincidence-cluster with updated Reco Pulses
      CrvCoincidenceCluster scrvcc(crvcc);
      scrvcc.SetCrvRecoPulses(pulses);
      scrvccc->push_back(scrvcc);
    }
    // Fill CrvIndex map
    uint16_t crvcount(0);
    for(auto crvindex : crvindices){
      crvdim->addElement(crvindex,crvcount++);
      // deep-copy the selected CrvDigis
      scrvdc->push_back(crvdc.at(crvindex));
    }
    // update digi indices in the pulses
    for(auto& crvrp : *scrvrpc) {
      auto& indices = crvrp.GetWaveformIndices();
      for(size_t iindex=0; iindex < indices.size(); ++iindex){
        indices[iindex] = crvdim->getCondensedIndex(indices[iindex]);
      }
    }
    // update reco count
    nrec._ncrvdigi = crvdc.size();
    // put new data in event
    event.put(std::move(scrvdc));
    event.put(std::move(scrvrpc));
    event.put(std::move(scrvccc));
    event.put(std::move(crvdim),"CrvDigiMap");
  }
  void SelectReco::fillCalo(art::Event& event, std::set<art::Ptr<CaloCluster> >& ccptrs, RecoCount& nrec) {
    auto cdch = event.getValidHandle<CaloDigiCollection>(cdct_);
    auto const& cdc = *cdch;
    auto ccch = event.getValidHandle<CaloClusterCollection>(ccct_);
    auto const& ccc = *ccch;

    std::unique_ptr<CaloDigiCollection> scdc(new CaloDigiCollection);
    // reco count
    nrec._ncalodigi = cdc.size();
    nrec._ncc = ccc.size();
    nrec._cce = 0.0;
    // loop over all the CaloClusters and mark the ones that are above energy for saving by adding their Ptrs to the list
    for(unsigned icc=0;icc < ccc.size(); icc++){
      auto const& cc = ccc[icc];
      nrec._cce += cc.energyDep();
      if(cc.energyDep() > ccme_){
        auto ccp = art::Ptr<CaloCluster>(ccch,icc);
        ccptrs.insert(ccp);
      }
    }
    // deep-copy CaloDigis from selected clusters
    for(auto const& ccptr : ccptrs) {
      for(auto const& cchptr : ccptr->caloHitsPtrVector()){
        for (auto const& rcdptr : cchptr->recoCaloDigis()){
          // deep-copy CaloDigis used in clusters
          scdc->push_back(*rcdptr->caloDigiPtr());
        }
      }
    }
    // put new data into event
    event.put(std::move(scdc));
  }

}
DEFINE_ART_MODULE(mu2e::SelectReco)
