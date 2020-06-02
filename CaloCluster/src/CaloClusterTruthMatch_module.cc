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
#include "MCDataProducts/inc/CaloDigiMC.hh"
#include "MCDataProducts/inc/CaloClusterMC.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

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
           caloDigiMCTruthToken_{consumes<CaloDigiMCTruthAssn>  (config().caloDigiMCTruthAssn())},
           diagLevel_           (config().diagLevel())
        {
           produces<CaloClusterMCCollection>();    
           produces<CaloClusterMCTruthAssn>();    
        }

	void produce(art::Event& e) override;


     private:
	void makeTruthMatch(art::Event&, CaloClusterMCCollection&,CaloClusterMCTruthAssn&);
        
        const art::ProductToken<CaloClusterCollection>  caloClusterToken_;
        const art::ProductToken<CaloDigiMCTruthAssn>    caloDigiMCTruthToken_;
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
           const auto& hits = cluster.caloCrystalHitsPtrVector();
           
           if (diagLevel_ > 1) std::cout<<"[CaloClusterTruthMatch] Inspect cluster diskId/energy/time "<<cluster.diskId()<<" "<<cluster.energyDep()<<" "<<cluster.time()<<std::endl;
           std::vector<CaloEDepMC> edeps;  
                     
           for (auto i=caloHitTruth.begin(), ie = caloHitTruth.end(); i !=ie; ++i)
           {	       
	        if (std::find(hits.begin(),hits.end(),i->first) == hits.end()) continue;
                const auto& digiMC = i->second;
                
                if (diagLevel_ > 1) std::cout<<"[CaloClusterTruthMatch] found hit in map "<<digiMC->nParticles()<<" "<<digiMC->time()<<std::endl;
                
                for (const auto& digiEdep : digiMC->energyDeposits())
                {
                    auto it = edeps.begin();
                    while (it != edeps.end()) {if (it->sim() == digiEdep.sim()) break;  ++it;}
                    
                    if (it!= edeps.end()) 
                    {
                        it->addEDep(digiEdep.energyDep());
                        it->addEDepG4(digiEdep.energyDepG4());
                        it->addTime(digiEdep.time());
                        it->addMom(digiEdep.momentumIn());
                    }
                    else 
                    {
                        edeps.emplace_back(CaloEDepMC(digiEdep.sim(),digiEdep.energyDep(),digiEdep.energyDepG4(),
                                                      digiEdep.time(),digiEdep.momentumIn(),digiEdep.rel()));                
                    }
                } 
                
           } 
           if (edeps.empty()) continue;
           
           std::sort(edeps.begin(),edeps.end(),[](const auto& a, const auto& b){return a.energyDep() > b.energyDep();});

           CaloClusterMCs.emplace_back(CaloClusterMC(std::move(edeps)));

           art::Ptr<CaloClusterMC> clusterMCPtr = art::Ptr<CaloClusterMC>(clusterMCProductID, CaloClusterMCs.size()-1, clusterMCProductGetter);             
           caloClusterTruthMatch.addSingle(clusterPtr,clusterMCPtr);
           
           totalEnergyMatched += clusterPtr->energyDep();++nMatched;
      }            	      
      
      if (diagLevel_ > 0) std::cout<<"[CaloClusterTruthMatch]  total clusters / energy matched = "<<nMatched<<" / "<<totalEnergyMatched<<std::endl;
 } 


 


}

using mu2e::CaloClusterTruthMatch;
DEFINE_ART_MODULE(CaloClusterTruthMatch);



