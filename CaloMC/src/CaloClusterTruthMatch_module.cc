//
// An EDProducer Module to match calo clusters to MC info
//
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"

#include "Offline/MCDataProducts/inc/CaloEDepMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include "Offline/MCDataProducts/inc/CaloClusterMC.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"

#include <iostream>
#include <string>
#include <map>
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
        void makeTruthMatch(art::Event&, CaloClusterMCCollection&,CaloClusterMCTruthAssn&);

        const art::ProductToken<CaloClusterCollection>  caloClusterToken_;
        const art::ProductToken<CaloHitMCTruthAssn>     caloHitMCTruthToken_;
        int                                             diagLevel_;
  };





  //--------------------------------------------------------------------
  void CaloClusterTruthMatch::produce(art::Event& event)
  {
      std::unique_ptr<CaloClusterMCCollection> caloClusterMCs(new CaloClusterMCCollection);
      std::unique_ptr<CaloClusterMCTruthAssn>  CaloClusterMCTruth(new CaloClusterMCTruthAssn);

      makeTruthMatch(event, *caloClusterMCs, *CaloClusterMCTruth);

      event.put(std::move(CaloClusterMCTruth));
      event.put(std::move(caloClusterMCs));
  }


  //--------------------------------------------------------------------
  void CaloClusterTruthMatch::makeTruthMatch(art::Event& event, CaloClusterMCCollection& caloClusterMCs,
                                             CaloClusterMCTruthAssn& caloClusterTruthMatch)
  {

      art::ProductID clusterMCProductID(event.getProductID<CaloClusterMCCollection>());
      const art::EDProductGetter* clusterMCProductGetter = event.productGetter(clusterMCProductID);

      const auto  caloClusterHandle = event.getValidHandle(caloClusterToken_);
      const auto& caloClusters(*caloClusterHandle);
      const auto* caloClusterBase = caloClusters.data();

      const auto  CaloHitMCHandle = event.getValidHandle(caloHitMCTruthToken_);
      const auto& caloHitTruth(*CaloHitMCHandle);

      double totalEnergyMatched(0);
      int nMatched(0);


      for (const auto& cluster : caloClusters)
      {
              const CaloCluster* thisCaloCluster = &cluster;
           size_t idx = (thisCaloCluster - caloClusterBase);
           art::Ptr<CaloCluster> clusterPtr = art::Ptr<CaloCluster>(caloClusterHandle,idx);
           const auto& hits = cluster.caloHitsPtrVector();

           if (diagLevel_ > 1) std::cout<<"[CaloClusterTruthMatch] Inspect cluster diskId/energy/time "<<cluster.diskID()<<" "<<cluster.energyDep()<<" "<<cluster.time()<<std::endl;
           std::vector<art::Ptr<CaloHitMC>> digis;

           for (auto i=caloHitTruth.begin(), ie = caloHitTruth.end(); i !=ie; ++i)
           {
                if (std::find(hits.begin(),hits.end(),i->first) == hits.end()) continue;
                const auto& digiMC = i->second;
                digis.push_back(digiMC);

                if (diagLevel_ > 1 && digiMC->nParticles()>0) std::cout<<"[CaloClusterTruthMatch] found hit in map "<<digiMC->nParticles()<<" "<<digiMC->time()<<std::endl;
           }

           std::sort(digis.begin(),digis.end(),[](const auto& a, const auto& b){return a->totalEnergyDep() > b->totalEnergyDep();});
           caloClusterMCs.emplace_back(CaloClusterMC(std::move(digis)));

           art::Ptr<CaloClusterMC> clusterMCPtr = art::Ptr<CaloClusterMC>(clusterMCProductID, caloClusterMCs.size()-1, clusterMCProductGetter);
           caloClusterTruthMatch.addSingle(clusterPtr,clusterMCPtr);

           totalEnergyMatched += clusterPtr->energyDep();
           ++nMatched;
      }

      if (diagLevel_ > 0) std::cout<<"[CaloClusterTruthMatch]  total clusters / energy matched = "<<nMatched<<" / "<<totalEnergyMatched<<std::endl;
 }

}

using mu2e::CaloClusterTruthMatch;
DEFINE_ART_MODULE(CaloClusterTruthMatch)



