//
// An EDProducer Module to match calo clusters to MC info
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/MCDataProducts/inc/CaloEDepMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include "Offline/MCDataProducts/inc/CaloClusterMC.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <vector>


namespace mu2e {

  class CaloClusterTruthMatch : public art::EDProducer {

     public:
         struct Config
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             fhicl::Atom<art::InputTag>  caloClusterCollection { Name("caloClusterCollection"), Comment("Name of calo cluster Collection") };
             fhicl::Atom<art::InputTag>  caloHitMCTruthAssn    { Name("caloHitMCTruthAssn"),    Comment("Name of caloHit - CaloHitMC Assn") };
             fhicl::Atom<int>            diagLevel             { Name("diagLevel"),             Comment("Diag Level"),0 };
         };

        explicit CaloClusterTruthMatch(const art::EDProducer::Table<Config>& config) :
           EDProducer{config},
           caloClusterToken_    {consumes<CaloClusterCollection>(config().caloClusterCollection())},
           caloHitMCTruthToken_ {consumes<CaloHitMCTruthAssn>   (config().caloHitMCTruthAssn())},
           diagLevel_           (config().diagLevel())
        {
           produces<CaloClusterMCCollection>();
           produces<CaloClusterMCTruthAssn>();
        }

        void produce(art::Event& e) override;


     private:
        void makeTruthMatch(art::Event&, CaloClusterMCCollection&, CaloClusterMCTruthAssn&) const;

        const art::ProductToken<CaloClusterCollection>  caloClusterToken_;
        const art::ProductToken<CaloHitMCTruthAssn>     caloHitMCTruthToken_;
        const int                                       diagLevel_;
  };



  //--------------------------------------------------------------------
  void CaloClusterTruthMatch::produce(art::Event& event)
  {
      auto caloClusterMCs   = std::make_unique<CaloClusterMCCollection>();
      auto caloClusterMCTruth = std::make_unique<CaloClusterMCTruthAssn>();

      makeTruthMatch(event, *caloClusterMCs, *caloClusterMCTruth);

      event.put(std::move(caloClusterMCTruth));
      event.put(std::move(caloClusterMCs));
  }


  //--------------------------------------------------------------------
  void CaloClusterTruthMatch::makeTruthMatch(art::Event& event, CaloClusterMCCollection& caloClusterMCs,
                                             CaloClusterMCTruthAssn& caloClusterTruthMatch) const
  {
      const art::ProductID        clusterMCProductID(event.getProductID<CaloClusterMCCollection>());
      const art::EDProductGetter* clusterMCProductGetter = event.productGetter(clusterMCProductID);

      const auto  caloClusterHandle = event.getValidHandle(caloClusterToken_);
      const auto& caloClusters(*caloClusterHandle);
      const auto& caloHitTruth = *event.getValidHandle(caloHitMCTruthToken_);

      // build the hit -> hitMC lookup once, instead of re-scanning the whole association per cluster
      std::map<art::Ptr<CaloHit>, art::Ptr<CaloHitMC>> hitToMC;
      for (const auto& assn : caloHitTruth) hitToMC.emplace(assn.first, assn.second);

      double totalEnergyMatched(0);
      int    nMatched(0);

      caloClusterMCs.reserve(caloClusters.size());

      for (std::size_t idx=0; idx<caloClusters.size(); ++idx) {
          const CaloCluster& cluster    = caloClusters[idx];
          const auto         clusterPtr = art::Ptr<CaloCluster>(caloClusterHandle,idx);

          if (diagLevel_ > 1) std::cout<<"[CaloClusterTruthMatch] Inspect cluster diskId/energy/time "
                                       <<cluster.diskID()<<" "<<cluster.energyDep()<<" "<<cluster.time()<<std::endl;

          // gather the CaloHitMC of every hit in this cluster that carries MC truth
          std::vector<art::Ptr<CaloHitMC>> digis;
          for (const auto& hitPtr : cluster.caloHitsPtrVector()) {
              const auto it = hitToMC.find(hitPtr);
              if (it == hitToMC.end()) continue;

              digis.push_back(it->second);
              if (diagLevel_ > 1 && it->second->nParticles()>0)
                 std::cout<<"[CaloClusterTruthMatch] found hit in map "<<it->second->nParticles()<<" "<<it->second->time()<<std::endl;
          }

          std::sort(digis.begin(),digis.end(),[](const auto& a, const auto& b){return a->totalEnergyDep() > b->totalEnergyDep();});
          caloClusterMCs.emplace_back(std::move(digis));

          const auto clusterMCPtr = art::Ptr<CaloClusterMC>(clusterMCProductID, caloClusterMCs.size()-1, clusterMCProductGetter);
          caloClusterTruthMatch.addSingle(clusterPtr,clusterMCPtr);

          totalEnergyMatched += cluster.energyDep();
          ++nMatched;
      }

      if (diagLevel_ > 0) std::cout<<"[CaloClusterTruthMatch]  total clusters / energy matched = "<<nMatched<<" / "<<totalEnergyMatched<<std::endl;
  }

}

DEFINE_ART_MODULE(mu2e::CaloClusterTruthMatch)
