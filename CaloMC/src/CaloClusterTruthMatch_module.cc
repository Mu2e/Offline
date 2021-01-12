//
// An EDProducer Module to match calo clusters to MC info
//
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"

#include "MCDataProducts/inc/CaloEDepMC.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "MCDataProducts/inc/CaloHitMC.hh"
#include "MCDataProducts/inc/CaloClusterMC.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"

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
             fhicl::Atom<art::InputTag>  caloDigiMCTruthAssn   { Name("caloDigiMCTruthAssn"),   Comment("Name of caloHit - caloDigiMC Assn") };
             fhicl::Atom<int>            diagLevel             { Name("diagLevel"),             Comment("Diag Level"),0 };
         };

        explicit CaloClusterTruthMatch(const art::EDProducer::Table<Config>& config) :
           EDProducer{config},
           caloClusterToken_    {consumes<CaloClusterCollection>(config().caloClusterCollection())},
           caloDigiMCTruthToken_{consumes<CaloHitMCTruthAssn>  (config().caloDigiMCTruthAssn())},
           diagLevel_           (config().diagLevel())
        {
           produces<CaloClusterMCCollection>();    
           produces<CaloClusterMCTruthAssn>();    
        }

	void produce(art::Event& e) override;


     private:
	void makeTruthMatch(art::Event&, CaloClusterMCCollection&,CaloClusterMCTruthAssn&);
        
        const art::ProductToken<CaloClusterCollection>  caloClusterToken_;
        const art::ProductToken<CaloHitMCTruthAssn>    caloDigiMCTruthToken_;
	int                                             diagLevel_;
  };





  //--------------------------------------------------------------------
  void CaloClusterTruthMatch::produce(art::Event& event)
  {
      std::unique_ptr<CaloClusterMCCollection> CaloClusterMCs(new CaloClusterMCCollection);
      std::unique_ptr<CaloClusterMCTruthAssn>  CaloClusterMCTruth(new CaloClusterMCTruthAssn);
   
      makeTruthMatch(event, *CaloClusterMCs, *CaloClusterMCTruth);

      event.put(std::move(CaloClusterMCTruth));
      event.put(std::move(CaloClusterMCs));
  } 

  
  //--------------------------------------------------------------------
  void CaloClusterTruthMatch::makeTruthMatch(art::Event& event, CaloClusterMCCollection& CaloClusterMCs, 
                                             CaloClusterMCTruthAssn& caloClusterTruthMatch)
  {
        
      art::ProductID clusterMCProductID(event.getProductID<CaloClusterMCCollection>());
      const art::EDProductGetter* clusterMCProductGetter = event.productGetter(clusterMCProductID);

      const auto  caloClusterHandle = event.getValidHandle(caloClusterToken_);
      const auto& caloClusters(*caloClusterHandle);
      const auto* caloClusterBase = &caloClusters.front();

      const auto  caloDigiMCHandle = event.getValidHandle(caloDigiMCTruthToken_);
      const auto& caloHitTruth(*caloDigiMCHandle);

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

                if (diagLevel_ > 1) std::cout<<"[CaloClusterTruthMatch] found hit in map "<<digiMC->nParticles()<<" "<<digiMC->time()<<std::endl;
           } 
           if (digis.empty()) continue;
           
           std::sort(digis.begin(),digis.end(),[](const auto& a, const auto& b){return a->totalEnergyDep() > b->totalEnergyDep();});
           CaloClusterMCs.emplace_back(CaloClusterMC(std::move(digis)));

           art::Ptr<CaloClusterMC> clusterMCPtr = art::Ptr<CaloClusterMC>(clusterMCProductID, CaloClusterMCs.size()-1, clusterMCProductGetter);             
           caloClusterTruthMatch.addSingle(clusterPtr,clusterMCPtr);
           
           totalEnergyMatched += clusterPtr->energyDep();
           ++nMatched;
      }            	      
      
      if (diagLevel_ > 0) std::cout<<"[CaloClusterTruthMatch]  total clusters / energy matched = "<<nMatched<<" / "<<totalEnergyMatched<<std::endl;
 } 

}

using mu2e::CaloClusterTruthMatch;
DEFINE_ART_MODULE(CaloClusterTruthMatch);



