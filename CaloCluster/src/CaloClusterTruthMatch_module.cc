//
// An EDProducer Module that matched caloCrystalHits to MC info
//
// Original author B. Echenard
//
// We match the CaloCluster to the caloShower (both time are folded), then extract the Simparticles
//

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"


// Mu2e includes.
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerCollection.hh"
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
	
	  _caloClusterModuleLabel  (pset.get<std::string>("caloClusterModuleLabel")), 
          _caloHitTruthModuleLabel (pset.get<std::string>("caloHitTruthModuleLabel")),
	  _diagLevel               (pset.get<int>        ("diagLevel",0))		  
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

	std::string  _caloClusterModuleLabel;   
        std::string  _caloHitTruthModuleLabel;
	int          _diagLevel;

        TH1F*  _hGenId;
        TH1F*  _hEner0;
        TH1F*  _hEner;
	
	
	void makeTruthMatch(CaloClusterMCTruthAssns &caloClusterTruthMatch, 
                            const CaloHitMCTruthAssns& caloHitTruth, 
			    const art::Handle<CaloClusterCollection> &CaloClusterHandle);

  };



  //--------------------------------------------------------------------
  void CaloClusterTruthMatch::beginJob()
  {
       art::ServiceHandle<art::TFileService> tfs;
       _hGenId    = tfs->make<TH1F>("hSimId",    "Sim gen Id",            150,    -10,  140);
       _hEner0    = tfs->make<TH1F>("hEner0",    "Signal cluster energy", 150,      0,  150);
       _hEner     = tfs->make<TH1F>("hEner",     "Signal cluster energy", 150,      0,  150);
  }




  //--------------------------------------------------------------------
  void CaloClusterTruthMatch::produce(art::Event& event)
  {
   
      art::Handle<CaloClusterCollection> caloClusterHandle;
      event.getByLabel(_caloClusterModuleLabel, caloClusterHandle);
 
      art::Handle<CaloHitMCTruthAssns> caloHitTruthHandle;
      event.getByLabel(_caloHitTruthModuleLabel, caloHitTruthHandle);
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
	    const auto& caloShowerPtr = caloHitTruth.data(i);
	    caloClusterTruth.addSingle(clusterPtr, sim, caloShowerPtr);

	    /*
	    if (_diagLevel > 0) 
	    {
	        _hEner0->Fill(thisCaloCluster->energyDep());
	       if (sim->genParticle()) _hGenId->Fill(sim->genParticle()->generatorId().id());
	       if (sim->genParticle() && sim->genParticle()->generatorId().id()==2) _hEner->Fill(thisCaloCluster->energyDep());
	       
	    }
	    */
            
	    if (_diagLevel > 2) std::cout<<"[CaloClusterTruthMatch]  matched crystal  id/time/Edep= "<<caloShowerPtr->crystalId()<<" / "<<caloShowerPtr->time()<<" / "<<caloShowerPtr->energy()
		                          <<"\t    hit in cluster id/time/Edep= "<<i->first->id()<<" / "<<i->first->time()<<" / "<<i->first->energyDep()<<std::endl;
       }
       
       
       if (_diagLevel > 0)
       {
           for (auto i=caloClusterTruth.begin(), ie = caloClusterTruth.end(); i !=ie; ++i)
	   {
	         const auto& cluster = i->first;
		 const auto& sim = i->second;
		 
	        _hEner0->Fill(cluster->energyDep());
	        if (sim->genParticle()) _hGenId->Fill(sim->genParticle()->generatorId().id());
	        if (sim->genParticle() && sim->genParticle()->generatorId().id()==2) _hEner->Fill(cluster->energyDep());
	   }
       }
                      
            	
  } 









 


}

using mu2e::CaloClusterTruthMatch;
DEFINE_ART_MODULE(CaloClusterTruthMatch);



