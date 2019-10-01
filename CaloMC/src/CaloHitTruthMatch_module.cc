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

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerSimCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthAssn.hh"
#include "MCDataProducts/inc/SimParticle.hh"

#include "TH2F.h"
#include "TFile.h"

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>


namespace mu2e {


  class CaloHitTruthMatch : public art::EDProducer {
  
     public:

        explicit CaloHitTruthMatch(fhicl::ParameterSet const& pset) :
          art::EDProducer{pset},
          // Parameters
          caloShowerSimModuleLabel_  (pset.get<std::string>("caloShowerSimModuleLabel")), 
          caloCrystalHitModuleLabel_ (pset.get<std::string>("caloCrystalHitModuleLabel")), 
          deltaTimeMinus_            (pset.get<double>     ("deltaTimeMinus")),  
          deltaTimePlus_             (pset.get<double>     ("deltaTimePlus")),   
          diagLevel_                 (pset.get<int>        ("diagLevel",0))                  
        {  

          produces<CaloHitMCTruthAssns>();    

        }

        virtual ~CaloHitTruthMatch() { }
        virtual void beginJob();

        void produce(art::Event& e);




     private:

        std::string       caloShowerSimModuleLabel_;   
        std::string       caloCrystalHitModuleLabel_;   
        double            deltaTimeMinus_;
        double            deltaTimePlus_;
        int               diagLevel_;


        //some diagnostic histograms
        TH1F*  hTime_;
        TH1F*  hTime2_;
        TH2F*  hTime2d_;
        TH2F*  hEner2d_;
        TH2F*  hEnerTime_;
        TH2F*  hdEdT_;
        TH1F*  hChi2_;


        void makeTruthMatch(const art::Handle<CaloShowerSimCollection> &caloShowerSimHandle, 
                            const art::Handle<CaloCrystalHitCollection> &CaloCrystalHitHandle,
                            CaloHitMCTruthAssns &caloTruthMatch);

  };







  //--------------------------------------------------------------------
  void CaloHitTruthMatch::beginJob()
  {
       if ( diagLevel_ > 2)
       {
	  art::ServiceHandle<art::TFileService> tfs;
	  hTime_     = tfs->make<TH1F>("hTime",    "delta Time",   2000, -20., 180);
	  hTime2_    = tfs->make<TH1F>("hTime2",   "delta Time",   2000, -20., 180);
	  hTime2d_   = tfs->make<TH2F>("hTime2d",  "Reco vs Gen time",  200,500,1700, 200,500,1700);
	  hEner2d_   = tfs->make<TH2F>("hEner2d",  "Reco vs gen Ener",  200,0,40,  200,0.,40);
	  hEnerTime_ = tfs->make<TH2F>("hTimeEner","delta Time vs Ener",500,0,100, 300,-50.,250);
	  hdEdT_     = tfs->make<TH2F>("hdEdt",    "delta Time vs delta Ener",100,-10,70, 170,-10.,160);
	  hChi2_     = tfs->make<TH1F>("hChi2",    "chi2 large dE",     50, 0., 10);
       }
  }



  //--------------------------------------------------------------------
  void CaloHitTruthMatch::produce(art::Event& event)
  {
   
      art::Handle<CaloShowerSimCollection> caloShowerSimHandle;
      event.getByLabel(caloShowerSimModuleLabel_, caloShowerSimHandle);

      art::Handle<CaloCrystalHitCollection> CaloCrystalHitHandle;
      event.getByLabel(caloCrystalHitModuleLabel_, CaloCrystalHitHandle);

      std::unique_ptr<CaloHitMCTruthAssns> caloHitMCTruth(new CaloHitMCTruthAssns);


      makeTruthMatch(caloShowerSimHandle,CaloCrystalHitHandle,*caloHitMCTruth);


      event.put(std::move(caloHitMCTruth));
  } 

  
  
  // first, sort the CaloShowers/ CaloHits into vector for each crystalId
  // next, for each crystal, sort the corresponding vectors by time;
  // finally, perform the association with the following rules:
  //   MCtime must be inside the window [recoTime-deltaTimeMinus, recoTime+deltaTimePlus] to be associated to RecoHit
  //   MCtime must not already be in the window of the next hit, in which case it is associate to this one
  //
  //--------------------------------------------------------------------
  void CaloHitTruthMatch::makeTruthMatch(const art::Handle<CaloShowerSimCollection> &caloShowerSimHandle, 
                                             const art::Handle<CaloCrystalHitCollection> &caloCrystalHitHandle,
                                             CaloHitMCTruthAssns &caloTruthMatch)
  {
        
      const CaloShowerSimCollection&  caloShowerSims(*caloShowerSimHandle);
      const CaloCrystalHitCollection& caloCrystalHits(*caloCrystalHitHandle);
 
      const CaloCrystalHit* caloCrystalHitBase = &caloCrystalHits.front();
      const CaloShowerSim*  caloShowerSimBase  = &caloShowerSims.front();

           
      std::map<int, std::vector<const CaloShowerSim*>> caloShowerSimsMap;
      for (auto const& caloShowerSim: caloShowerSims) caloShowerSimsMap[caloShowerSim.crystalId()].push_back(&caloShowerSim);
      
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
          
          std::vector<const CaloShowerSim*>& caloShowerSims = caloShowerSimsMap[crystalId];
          std::sort(caloShowerSims.begin(),caloShowerSims.end(), [](auto const a, auto const b){return a->time() < b->time();});
          
          auto showerIt    = caloShowerSims.begin();
          auto showerItEnd = caloShowerSims.end();
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
                     hTime2d_->Fill((*showerIt)->time(),(*hitIt)->time());
                     hEner2d_->Fill((*showerIt)->energy(),(*hitIt)->energyDep());
		     hdEdT_->Fill((*hitIt)->energyDep()-(*showerIt)->energy(),(*showerIt)->time()-(*hitIt)->time());
		    
		     if ((*showerIt)->energy() > 5) hTime2_->Fill((*showerIt)->time()-(*hitIt)->time()); 
		    
		     double deltaE = std::abs((*showerIt)->energy()-(*hitIt)->energyDep());
		     if (deltaE > 5 && (*showerIt)->energy() > 5)
		     hChi2_->Fill((*hitIt)->recoCaloDigis().at(0)->chi2()/(*hitIt)->recoCaloDigis().at(0)->ndf());		 
		 }
                 
                 // add shower to the match
                 const CaloCrystalHit* hit      = *hitIt;
                 const CaloShowerSim* showerSim = *showerIt;
                 size_t idxHit                  = (hit - caloCrystalHitBase);
                 size_t idxShower               = (showerSim - caloShowerSimBase);               
                 auto hitPtr                    = art::Ptr<CaloCrystalHit>(caloCrystalHitHandle,idxHit);               
                 auto  ShowerSimPtr             = art::Ptr<CaloShowerSim>(caloShowerSimHandle,idxShower);
                 
		 caloTruthMatch.addSingle(hitPtr, showerSim->sim(), ShowerSimPtr);

                 if (diagLevel_ > 3) std::cout<<"[CaloHitTruthMatch]  matched shower id/time/energyDep()= "<<showerSim->crystalId()
                                              <<" / "<<showerSim->time()<<" / "<<showerSim->energy()
                                              <<"\t    hit  id/time/energyDep()= "<<hit->id()<<" / "<<hit->time()<<" / "<<hit->energyDep()<<std::endl;

                 ++showerIt;                
             }
             if (hitIsMatched) {totalEnergyMatched += (*hitIt)->energyDep();++nMatched;}

             ++hitIt;
          }      
      }      
      
      if (diagLevel_ > 0) std::cout<<"[CaloHitTruthMatch]  total particles / energy matched = "<<nMatched<<" / "<<totalEnergyMatched<<std::endl;

  } 




}

using mu2e::CaloHitTruthMatch;
DEFINE_ART_MODULE(CaloHitTruthMatch);



