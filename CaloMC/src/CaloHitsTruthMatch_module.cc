//
// An EDProducer Module that matched caloCrystalHits to MC info
//
// Original author B. Echenard
//
// A CaloCrystalHit is made of one or more CaloRecoDigis, each of which has been extracted from 
// the caloDigi, themselves created from the caloShowers.
// A caloShower is made of one or more caloShowerStepMC (including time and energy smearing), and each 
// caloShowerStepMC corresponds to one SimParticle 
//
// So we match the CaloCrystalHit to the caloShower (both time are folded), then extract the generated info 
// from caloShowerStepMC and SimParticles


// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>

// ROOT includes
#include "TH2F.h"
#include "TFile.h"


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
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthAssn.hh"
#include "MCDataProducts/inc/SimParticle.hh"




namespace mu2e {


  class CaloHitsTruthMatch : public art::EDProducer {
  
     public:

        explicit CaloHitsTruthMatch(fhicl::ParameterSet const& pset) :

          // Parameters
          caloShowerModuleLabel_     (pset.get<std::string>("caloShowerModuleLabel")), 
          caloCrystalHitModuleLabel_ (pset.get<std::string>("caloCrystalHitModuleLabel")), 
          deltaTimeMinus_            (pset.get<double>     ("deltaTimeMinus")),  
          deltaTimePlus_             (pset.get<double>     ("deltaTimePlus")),   
          diagLevel_                 (pset.get<int>        ("diagLevel",0))                  
        {  

          produces<CaloHitMCTruthAssns>();    

        }

        virtual ~CaloHitsTruthMatch() { }
        virtual void beginJob();

        void produce(art::Event& e);




     private:

        std::string       caloShowerModuleLabel_;   
        std::string       caloCrystalHitModuleLabel_;   
        double            deltaTimeMinus_;
        double            deltaTimePlus_;
        int               diagLevel_;


        //some diagnostic histograms
        TH1F*  hTime_;
        TH1F*  hTime2_;
        TH2F*  hEnerTime_;


        void makeTruthMatch(const art::Handle<CaloShowerCollection> &caloShowerHandle, 
                            const art::Handle<CaloCrystalHitCollection> &CaloCrystalHitHandle,
                            CaloHitMCTruthAssns &caloTruthMatch);

  };







  //--------------------------------------------------------------------
  void CaloHitsTruthMatch::beginJob()
  {
       art::ServiceHandle<art::TFileService> tfs;
       hTime_     = tfs->make<TH1F>("hTime",    "delta Time",   2000, -20., 180);
       hTime2_    = tfs->make<TH1F>("hTime2",   "delta Time",   2000, -20., 180);
       hEnerTime_ = tfs->make<TH2F>("hTimeEner","delta Time vs Ener",500,0,100,300,-50.,250);
  }



  //--------------------------------------------------------------------
  void CaloHitsTruthMatch::produce(art::Event& event)
  {
   
      art::Handle<CaloShowerCollection> caloShowerHandle;
      event.getByLabel(caloShowerModuleLabel_, caloShowerHandle);

      art::Handle<CaloCrystalHitCollection> CaloCrystalHitHandle;
      event.getByLabel(caloCrystalHitModuleLabel_, CaloCrystalHitHandle);

      std::unique_ptr<CaloHitMCTruthAssns> caloHitMCTruth(new CaloHitMCTruthAssns);


      makeTruthMatch(caloShowerHandle,CaloCrystalHitHandle,*caloHitMCTruth);


      event.put(std::move(caloHitMCTruth));
  } 

  
  
  // first, sort the CaloShowers/ CaloHits into vector for each crystalId
  // next, for each crystal, sort the corresponding vectors by time;
  // finally, perform the association with the following rules:
  //   MCtime must be inside the window [recoTime-deltaTimeMinus, recoTime+deltaTimePlus] to be associated to RecoHit
  //   MCtime must not already be in the window of the next hit, in which case it is associate to this one
  //
  //--------------------------------------------------------------------
  void CaloHitsTruthMatch::makeTruthMatch(const art::Handle<CaloShowerCollection> &caloShowerHandle, 
                                             const art::Handle<CaloCrystalHitCollection> &caloCrystalHitHandle,
                                             CaloHitMCTruthAssns &caloTruthMatch)
  {
        
      const CaloShowerCollection&     caloShowers(*caloShowerHandle);
      const CaloCrystalHitCollection& caloCrystalHits(*caloCrystalHitHandle);
 
      const CaloCrystalHit* caloCrystalHitBase = &caloCrystalHits.front();
      const CaloShower*     caloShowerBase     = &caloShowers.front();

           
      std::map<int, std::vector<const CaloShower*>> caloShowersMap;
      for (auto const& caloShower: caloShowers) caloShowersMap[caloShower.crystalId()].push_back(&caloShower);
      
      std::map<int, std::vector<const CaloCrystalHit*>> caloHitMap;
      for (auto const& CaloCrystalHit: caloCrystalHits) caloHitMap[CaloCrystalHit.id()].push_back(&CaloCrystalHit);
      
      int nMatched(0);
      double totalEnergyMatched(0); 
      for (auto &kv : caloHitMap)
      {
          int crystalId = kv.first;

          std::vector<const CaloCrystalHit*> &caloHits = kv.second;          
          if (!caloHits.size()) continue;
          std::sort(caloHits.begin(),caloHits.end(), [](auto const a, auto const b){return a->time() < b->time();});
          
          std::vector<const CaloShower*>& caloShowers = caloShowersMap[crystalId];
          std::sort(caloShowers.begin(),caloShowers.end(), [](auto const a, auto const b){return a->time() < b->time();});
          
          auto showerIt    = caloShowers.begin();
          auto showerItEnd = caloShowers.end();
          auto hitIt       = caloHits.begin();
          auto hitItEnd    = caloHits.end();



          while( hitIt != hitItEnd)
          {            
             auto hitNextIt = std::next(hitIt);
             bool hitIsMatched(false);

             //forward until we reach the recoHit time;
             while (showerIt != showerItEnd && ( (*showerIt)->time() < (*hitIt)->time() - deltaTimeMinus_) ) ++showerIt; 
                          
             //loop as long as the shower time is witthin the recoTime window
             while(showerIt != showerItEnd && ( (*showerIt)->time() < (*hitIt)->time() + deltaTimePlus_) )
             {
                 //check if we're already inside the next hit time window
                 if (hitNextIt != hitItEnd && ( (*showerIt)->time() > (*hitNextIt)->time() - deltaTimeMinus_) ) break;
                 hitIsMatched = true;
                 
                 if (diagLevel_ > 2)
                 {
                    hTime_->Fill((*showerIt)->time()-(*hitIt)->time());
                    hEnerTime_->Fill((*showerIt)->energy(),(*showerIt)->time()-(*hitIt)->time());
		    if ((*showerIt)->energy() > 5) hTime2_->Fill((*showerIt)->time()-(*hitIt)->time()); 
		 }
                 
                 // add shower to the match
                 const CaloCrystalHit* hit = *hitIt;
                 const CaloShower* shower  = *showerIt;
                 size_t idxHit             = (hit - caloCrystalHitBase);
                 size_t idxShower          = (shower - caloShowerBase);               
                 auto hitPtr               = art::Ptr<CaloCrystalHit>(caloCrystalHitHandle,idxHit);               
                 auto  ShowerPtr           = art::Ptr<CaloShower>(caloShowerHandle,idxShower);

                 std::set<art::Ptr<SimParticle>> sims;
                 for (auto const& step : shower->caloShowerSteps()) sims.insert(step->simParticle());
                 for (auto &sim : sims) caloTruthMatch.addSingle(hitPtr, sim, ShowerPtr);

                 if (diagLevel_ > 3) std::cout<<"[CaloHitsTruthMatch]  matched shower id/time/energyDep()= "<<shower->crystalId()<<" / "<<shower->time()<<" / "<<shower->energy()
                                              <<"\t    hit  id/time/energyDep()= "<<hit->id()<<" / "<<hit->time()<<" / "<<hit->energyDep()<<std::endl;

                 ++showerIt;                
             }
             if (hitIsMatched) {totalEnergyMatched += (*hitIt)->energyDep();++nMatched;}

             ++hitIt;
          }      
      }      
      
      if (diagLevel_ > 0) std::cout<<"[CaloHitsTruthMatch]  total particles / energy matched = "<<nMatched<<" / "<<totalEnergyMatched<<std::endl;

  } 









 


}

using mu2e::CaloHitsTruthMatch;
DEFINE_ART_MODULE(CaloHitsTruthMatch);



