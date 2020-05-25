//
// An EDProducer Module that matched caloCrystalHits to MC info
//
// Original author B. Echenard
//
// We match the CaloCluster to the CaloShowerSim (both time are folded), then extract the Simparticles
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "MCDataProducts/inc/CaloEDepMC.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "MCDataProducts/inc/CaloDigiMCCollection.hh"
#include "MCDataProducts/inc/CaloClusterMCCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"

#include "TH2F.h"
#include "TFile.h"

#include <iostream>
#include <string>
#include <cmath>
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
           caloClusterToken_    {consumes<CaloClusterCollection> (config().caloClusterCollection())},
           caloDigiMCTruthToken_{consumes<CaloDigiMCTruthAssn>   (config().caloDigiMCTruthAssn())},
           diagLevel_           (config().diagLevel())
        {
           produces<CaloClusterNewMCCollection>();    
           produces<CaloClusterNewMCTruthAssn>();    
        }

	virtual ~CaloClusterTruthMatch() {}

	virtual void beginJob();
	void produce(art::Event& e);




     private:
        const art::ProductToken<CaloClusterCollection>  caloClusterToken_;
        const art::ProductToken<CaloDigiMCTruthAssn>    caloDigiMCTruthToken_;
	int                                             diagLevel_;

        TH1F*  hGenId_;
        TH1F*  hEner0_;
        TH1F*  hEner_;
	
	
	void makeTruthMatch(art::Event& event, CaloClusterNewMCCollection& CaloClusterNewMCs,CaloClusterNewMCTruthAssn& caloClusterTruthMatch);

  };



  //--------------------------------------------------------------------
  void CaloClusterTruthMatch::beginJob()
  {
      if (diagLevel_ > 2)
      {
          art::ServiceHandle<art::TFileService> tfs;
          hGenId_    = tfs->make<TH1F>("hSimId",    "Sim gen Id",            150,    -10,  140);
          hEner0_    = tfs->make<TH1F>("hEner0",    "Signal cluster energy", 150,      0,  150);
          hEner_     = tfs->make<TH1F>("hEner",     "Signal cluster energy", 150,      0,  150);
      }
  }




  //--------------------------------------------------------------------
  void CaloClusterTruthMatch::produce(art::Event& event)
  {
      std::unique_ptr<CaloClusterNewMCCollection> CaloClusterNewMCs(new CaloClusterNewMCCollection);
      std::unique_ptr<CaloClusterNewMCTruthAssn>  CaloClusterNewMCTruth(new CaloClusterNewMCTruthAssn);
   
      makeTruthMatch(event, *CaloClusterNewMCs, *CaloClusterNewMCTruth);

      event.put(std::move(CaloClusterNewMCTruth));
      event.put(std::move(CaloClusterNewMCs));
  } 

  
  //--------------------------------------------------------------------
  void CaloClusterTruthMatch::makeTruthMatch(art::Event& event, CaloClusterNewMCCollection& CaloClusterNewMCs, 
                                             CaloClusterNewMCTruthAssn& caloClusterTruthMatch)
  {
        
      art::ProductID clusterMCProductID(event.getProductID<CaloClusterNewMCCollection>());
      const art::EDProductGetter* clusterMCProductGetter = event.productGetter(clusterMCProductID);

      const auto  caloClusterHandle = event.getValidHandle(caloClusterToken_);
      const auto& caloClusters(*caloClusterHandle);
      const auto* caloClusterBase = &caloClusters.front();

      const auto  caloDigiMCHandle = event.getValidHandle(caloDigiMCTruthToken_);
      const auto& caloHitTruth(*caloDigiMCHandle);


      for (const auto& cluster : caloClusters)
      {
   	   const CaloCluster* thisCaloCluster = &cluster;
	   size_t idx = (thisCaloCluster - caloClusterBase);
	   art::Ptr<CaloCluster> clusterPtr = art::Ptr<CaloCluster>(caloClusterHandle,idx);           
           const auto& hits = cluster.caloCrystalHitsPtrVector();
           
           std::vector<CaloEDepMC> edeps;  
                     
           for (auto i=caloHitTruth.begin(), ie = caloHitTruth.end(); i !=ie; ++i)
           {	       
	        if (std::find(hits.begin(),hits.end(),i->first) == hits.end()) continue;
                const auto& digiMC = i->second;

                for (const auto& digiEdep : digiMC->energyDeposits())
                {
                    auto it = edeps.begin();
                    while (it != edeps.end()) {if (it->sim() == digiEdep.sim()) break;  ++it;}

                    if (it!= edeps.end()) 
                      it->set(it->eDep()+digiEdep.eDep(), it->eDepG4()+digiEdep.eDepG4(),
                              std::min(it->time(),digiEdep.time()), std::max(digiEdep.momentumIn(),it->momentumIn())); 
                    else 
                      edeps.emplace_back(CaloEDepMC(digiEdep.sim(),digiEdep.eDep(),digiEdep.eDepG4(),digiEdep.time(),digiEdep.momentumIn()));
                } 

                
           } 
           std::sort(edeps.begin(),edeps.end(),[](const auto& a, const auto& b){return a.eDep() > b.eDep();});

           CaloClusterNewMCs.emplace_back(CaloClusterNewMC(std::move(edeps)));

           art::Ptr<CaloClusterNewMC> clusterMCPtr = art::Ptr<CaloClusterNewMC>(clusterMCProductID, CaloClusterNewMCs.size(), clusterMCProductGetter);             
           caloClusterTruthMatch.addSingle(clusterPtr,clusterMCPtr);

      }

       /*
       if (diagLevel_ > 2)
       {
           for (auto i=caloClusterTruthMatch.begin(), ie = caloClusterTruthMatch.end(); i !=ie; ++i)
	   {
	       const auto& cluster       = i->first;
	       const auto& clusterDigiMC = i->second;
	       hEner0_->Fill(clusterDigiMC->energyDepTot());
	       //if (sim->genParticle()) hGenId_->Fill(clusterDigiMC->sim()->genParticle()->generatorId().id());
	       //if (sim->genParticle() && sim->genParticle()->generatorId().id()==2) hEner_->Fill(cluster->energyDep());
	   }
       }
       */
                      
            	
  } 









 


}

using mu2e::CaloClusterTruthMatch;
DEFINE_ART_MODULE(CaloClusterTruthMatch);



