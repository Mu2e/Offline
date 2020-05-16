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
#include "fhiclcpp/types/Atom.h"

#include "CaloMC/inc/CaloSimSummary.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerSimCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthAssn.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/CaloDigiMC.hh"

#include "TH2F.h"
#include "TFile.h"

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>



namespace mu2e {

  class CaloHitTruthMatch : public art::EDProducer 
  {
      public:
          struct Config 
          {
              using Name    = fhicl::Name;
              using Comment = fhicl::Comment;

              fhicl::Atom<art::InputTag>  caloShowerSimCollection  { Name("caloShowerSimCollection"),  Comment("Name of caloShowerSim Collection") };
              fhicl::Atom<art::InputTag>  caloCrystalHitCollection { Name("caloCrystalHitCollection"), Comment("Name of caloCrystalHit collection") };
              fhicl::Atom<double>         deltaTimeMinus           { Name("deltaTimeMinus"),           Comment("Maximum time before hit to include MC hit") };
              fhicl::Atom<double>         deltaTimePlus            { Name("deltaTimePlus"),            Comment("Maximum time after hit to include MC hit") };
              fhicl::Atom<int>            diagLevel                { Name("diagLevel"),                Comment("Diag Level"),0 };
          };


         explicit CaloHitTruthMatch(const art::EDProducer::Table<Config>& config) :
            EDProducer{config},
            caloShowerSimToken_ {consumes<CaloShowerSimCollection> (config().caloShowerSimCollection())},
            caloCrystalHitToken_{consumes<CaloCrystalHitCollection>(config().caloCrystalHitCollection())},
            deltaTimeMinus_     (config().deltaTimeMinus()),
            deltaTimePlus_      (config().deltaTimePlus()),
            diagLevel_          (config().diagLevel())
         {
             produces<CaloHitMCTruthAssns>();    
         }

         virtual void beginJob();
         void produce(art::Event& e);


      private:
          using SimParticlePtr = art::Ptr<SimParticle>;

          void makeTruthMatch(const art::ValidHandle<CaloShowerSimCollection>& caloShowerSimHandle, 
                              const art::ValidHandle<CaloCrystalHitCollection>& caloCrystalHitHandle,
                              CaloHitMCTruthAssns &caloTruthMatch);
          void diag(const CaloShowerSim* shower, const CaloCrystalHit* hit);


          const art::ProductToken<CaloShowerSimCollection>  caloShowerSimToken_;
          const art::ProductToken<CaloCrystalHitCollection> caloCrystalHitToken_;
          double  deltaTimeMinus_;
          double  deltaTimePlus_;
          int     diagLevel_;

          //some diagnostic histograms
          TH1F*  hTime_;
          TH1F*  hTime2_;
          TH2F*  hTime2d_;
          TH2F*  hEner2d_;
          TH2F*  hEnerTime_;
          TH2F*  hdEdT_;
          TH1F*  hChi2_;
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
      const auto caloShowerSimHandle  = event.getValidHandle(caloShowerSimToken_);
      const auto caloCrystalHitHandle = event.getValidHandle(caloCrystalHitToken_);
      std::unique_ptr<CaloHitMCTruthAssns> caloHitMCTruth(new CaloHitMCTruthAssns);

      makeTruthMatch(caloShowerSimHandle,caloCrystalHitHandle,*caloHitMCTruth);

      event.put(std::move(caloHitMCTruth));
  } 

  
  
  // first, sort the CaloShowers/ CaloHits into vector for each crystalId
  // next, for each crystal, sort the corresponding vectors by time;
  // finally, perform the association with the following rules:
  //   MCtime must be inside the window [recoTime-deltaTimeMinus, recoTime+deltaTimePlus] to be associated to RecoHit
  //   MCtime must not already be in the window of the next hit, in which case it is associate to this one
  //
  //--------------------------------------------------------------------
  void CaloHitTruthMatch::makeTruthMatch(const art::ValidHandle<CaloShowerSimCollection>& caloShowerSimHandle, 
                                         const art::ValidHandle<CaloCrystalHitCollection>& caloCrystalHitHandle,
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



          while (hitIt != hitItEnd)
          {            
             auto hitNextIt = std::next(hitIt);
             bool hitIsMatched(false);
             CaloSimSummary simMap;

             const CaloCrystalHit* hit = *hitIt;
             size_t idxHit             = (hit - caloCrystalHitBase);
             auto hitPtr               = art::Ptr<CaloCrystalHit>(caloCrystalHitHandle,idxHit);               
                         
             
             //forward until we reach the recoHit time;
             while (showerIt != showerItEnd && ( (*showerIt)->time() < (*hitIt)->time() - deltaTimeMinus_) ) ++showerIt; 
                          
             //loop as long as the shower time is witthin the recoTime window
             while(showerIt != showerItEnd && ( (*showerIt)->time() < (*hitIt)->time() + deltaTimePlus_) )
             {
                 //check if we're already inside the next hit time window
                 if (hitNextIt != hitItEnd && ( (*showerIt)->time() > (*hitNextIt)->time() - deltaTimeMinus_) ) break;
                 hitIsMatched = true;
                 
//Maybe out a check on energy to include the CaloShowerSim in the hit                 
                                  
                 const CaloShowerSim* showerSim(*showerIt);
                 size_t idxShower   = (showerSim - caloShowerSimBase);               
                 auto  ShowerSimPtr = art::Ptr<CaloShowerSim>(caloShowerSimHandle,idxShower);
                 float energyDep    = showerSim->energyDep();
                 float energyDepG4  = showerSim->energyDepG4();
                 float time         = showerSim->time();
                 float pIn          = showerSim->momentumIn();
                 
                 //add shower to SimParticle summary
                 simMap.add(showerSim->sim(),energyDep,energyDepG4,time,hit->time(),pIn);                 
                 
                 // add shower to the detailed truthMatch
                 caloTruthMatch.addSingle(hitPtr, showerSim->sim(), ShowerSimPtr);
                                 
                                                  
                 if (diagLevel_ > 2) diag(showerSim,hit);
                 if (diagLevel_ > 3) std::cout<<"[CaloHitTruthMatch]  matched shower id/time/energyDep()= "<<showerSim->crystalId()
                                              <<" / "<<showerSim->time()<<" / "<<showerSim->energyDep()
                                              <<"\t    hit  id/time/energyDep()= "<<hit->id()<<" / "<<hit->time()<<" / "<<hit->energyDep()<<std::endl;
                 ++showerIt;                
             }
             
             CaloDigiMC digiMC(simMap.sims(),simMap.energyDep(),simMap.energyDepG4(),simMap.time(),simMap.momentumIn(),simMap.isConversion());
             //need to push back in collection and add truth match             
             
             if (hitIsMatched) {totalEnergyMatched += (*hitIt)->energyDep();++nMatched;}
             ++hitIt;
          }      
      }      
      
      if (diagLevel_ > 0) std::cout<<"[CaloHitTruthMatch]  total particles / energy matched = "<<nMatched<<" / "<<totalEnergyMatched<<std::endl;




  } 


  void CaloHitTruthMatch::diag(const CaloShowerSim* shower, const CaloCrystalHit* hit)
  {
      hTime_->Fill(shower->time()-hit->time());
      hEnerTime_->Fill(shower->energyDep(),shower->time()-hit->time());
      hTime2d_->Fill(shower->time(),hit->time());
      hEner2d_->Fill(shower->energyDep(),hit->energyDep());
      hdEdT_->Fill(hit->energyDep()-shower->energyDep(),shower->time()-hit->time());

      if (shower->energyDep() > 5) hTime2_->Fill(shower->time()-hit->time()); 

      double deltaE = std::abs(shower->energyDep()-hit->energyDep());
      if (deltaE > 5 && shower->energyDep() > 5)
      hChi2_->Fill(hit->recoCaloDigis().at(0)->chi2()/hit->recoCaloDigis().at(0)->ndf());		   
  }


}

using mu2e::CaloHitTruthMatch;
DEFINE_ART_MODULE(CaloHitTruthMatch);



