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

#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerSimCollection.hh"
#include "MCDataProducts/inc/CaloClusterMCTruthAssn.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/CaloHitMCTruthAssn.hh"


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

	explicit CaloClusterTruthMatch(fhicl::ParameterSet const& pset) :
          art::EDProducer{pset},	
	  caloClusterModuleLabel_  (pset.get<std::string>("caloClusterModuleLabel")), 
          caloHitTruthModuleLabel_ (pset.get<std::string>("caloHitTruthModuleLabel")),
	  diagLevel_               (pset.get<int>        ("diagLevel",0))		  
	{  

	  produces<CaloClusterMCTruthAssns>();    

	}

	virtual ~CaloClusterTruthMatch() {}

	virtual void beginJob();
	void produce(art::Event& e);




     private:
        
	typedef art::Ptr<SimParticle>    SimParticlePtr;
	typedef art::Ptr<CaloCrystalHit> CaloCrystalHitPtr;
	typedef art::Ptr<CaloCluster>    CaloClusterPtr;

	std::string  caloClusterModuleLabel_;   
        std::string  caloHitTruthModuleLabel_;
	int          diagLevel_;

        TH1F*  hGenId_;
        TH1F*  hEner0_;
        TH1F*  hEner_;
	
	
	void makeTruthMatch(CaloClusterMCTruthAssns &caloClusterTruthMatch, 
                            const CaloHitMCTruthAssns& caloHitTruth, 
			    const art::Handle<CaloClusterCollection> &CaloClusterHandle);

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
   
      art::Handle<CaloClusterCollection> caloClusterHandle;
      event.getByLabel(caloClusterModuleLabel_, caloClusterHandle);
 
      art::Handle<CaloHitMCTruthAssns> caloHitTruthHandle;
      event.getByLabel(caloHitTruthModuleLabel_, caloHitTruthHandle);
      const CaloHitMCTruthAssns& caloHitTruth(*caloHitTruthHandle);

      std::unique_ptr<CaloClusterMCTruthAssns> caloClusterTruth(new CaloClusterMCTruthAssns);
 
 
      makeTruthMatch(*caloClusterTruth, caloHitTruth, caloClusterHandle);


      event.put(std::move(caloClusterTruth));
  } 

  
  
  //--------------------------------------------------------------------
  void CaloClusterTruthMatch::makeTruthMatch(CaloClusterMCTruthAssns& caloClusterTruth, 
                                                const CaloHitMCTruthAssns& caloHitTruth, 
						const art::Handle<CaloClusterCollection>& CaloClusterHandle)
  {
        
       const CaloClusterCollection& caloClusters(*CaloClusterHandle);
       const CaloCluster*           caloClusterBase = &caloClusters.front();


       std::map<CaloCrystalHitPtr, const CaloCluster*> hitToCluster;

       for (const auto& caloCluster : caloClusters)
	   for (const auto& caloCrystalHit : caloCluster.caloCrystalHitsPtrVector()) hitToCluster[caloCrystalHit] = &caloCluster;


       for (auto i=caloHitTruth.begin(), ie = caloHitTruth.end(); i !=ie; ++i)
       {	       
	    auto hitToClusterIt = hitToCluster.find(i->first);
	    if (hitToClusterIt == hitToCluster.end()) continue;

   	    const CaloCluster* thisCaloCluster = hitToClusterIt->second;
	    size_t idx = (thisCaloCluster - caloClusterBase);
	    art::Ptr<CaloCluster> clusterPtr = art::Ptr<CaloCluster>(CaloClusterHandle,idx);

	    const auto& sim = i->second;
	    const auto& CaloShowerSimPtr = caloHitTruth.data(i);
	    caloClusterTruth.addSingle(clusterPtr, sim, CaloShowerSimPtr);

            
	    if (diagLevel_ > 2) std::cout<<"[CaloClusterTruthMatch]  matched crystal  id/time/Edep= "<<CaloShowerSimPtr->crystalId()<<" / "<<CaloShowerSimPtr->time()<<" / "<<CaloShowerSimPtr->energy()
		                          <<"\t    hit in cluster id/time/Edep= "<<i->first->id()<<" / "<<i->first->time()<<" / "<<i->first->energyDep()<<std::endl;
       }
       
       
       if (diagLevel_ > 0)
       {
           for (auto i=caloClusterTruth.begin(), ie = caloClusterTruth.end(); i !=ie; ++i)
	   {
	         const auto& cluster = i->first;
		 const auto& sim = i->second;
		 
	         if (diagLevel_ > 2)
		 {
		    hEner0_->Fill(cluster->energyDep());
	            if (sim->genParticle()) hGenId_->Fill(sim->genParticle()->generatorId().id());
	            if (sim->genParticle() && sim->genParticle()->generatorId().id()==2) hEner_->Fill(cluster->energyDep());
	         }
	   }
       }
                      
            	
  } 









 


}

using mu2e::CaloClusterTruthMatch;
DEFINE_ART_MODULE(CaloClusterTruthMatch);



