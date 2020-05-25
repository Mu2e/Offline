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

#include "MCDataProducts/inc/CaloEDepMC.hh"
#include "MCDataProducts/inc/CaloDigiMCCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerSimCollection.hh"
#include "MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"

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
             fhicl::Atom<bool>           fillDetailedMC           { Name("fillDetailedMC"),           Comment("Fill SimParticle - SimShower Assn map")};
             fhicl::Atom<int>            diagLevel                { Name("diagLevel"),                Comment("Diag Level"),0 };
         };


        explicit CaloHitTruthMatch(const art::EDProducer::Table<Config>& config) :
           EDProducer{config},
           caloShowerSimToken_ {consumes<CaloShowerSimCollection> (config().caloShowerSimCollection())},
           caloCrystalHitToken_{consumes<CaloCrystalHitCollection>(config().caloCrystalHitCollection())},
           deltaTimeMinus_     (config().deltaTimeMinus()),
           deltaTimePlus_      (config().deltaTimePlus()),
           fillDetailedMC_     (config().fillDetailedMC()),
           diagLevel_          (config().diagLevel())
        {
            produces<CaloDigiMCCollection>();    
            produces<CaloDigiMCTruthAssn>();    
            if (fillDetailedMC_) produces<CaloShowerMCTruthAssn>();    
        }

        virtual void beginJob();
        void produce(art::Event& e);


      private:
         using SimParticlePtr = art::Ptr<SimParticle>;

         void makeTruthMatch(art::Event& event, CaloDigiMCCollection& caloDigiMCs,
                             CaloDigiMCTruthAssn& caloDigiTruthMatch, CaloShowerMCTruthAssn &caloShowerTruthMatch);
         void diag(const CaloShowerSim* shower, const CaloCrystalHit* hit);


         const art::ProductToken<CaloShowerSimCollection>  caloShowerSimToken_;
         const art::ProductToken<CaloCrystalHitCollection> caloCrystalHitToken_;
         double  deltaTimeMinus_;
         double  deltaTimePlus_;
         bool    fillDetailedMC_;
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
      std::unique_ptr<CaloDigiMCCollection>  caloDigiMCs(new CaloDigiMCCollection);
      std::unique_ptr<CaloDigiMCTruthAssn>   caloDigiMCTruth(new CaloDigiMCTruthAssn);
      std::unique_ptr<CaloShowerMCTruthAssn> caloHitMCTruth(new CaloShowerMCTruthAssn);

      makeTruthMatch(event, *caloDigiMCs, *caloDigiMCTruth, *caloHitMCTruth);

      event.put(std::move(caloDigiMCTruth));
      event.put(std::move(caloDigiMCs));
      if (fillDetailedMC_) event.put(std::move(caloHitMCTruth));
  } 

  
  
  // first, sort the CaloShowers/ CaloHits into vector for each crystalId
  // next, for each crystal, sort the corresponding vectors by time;
  // finally, perform the association with the following rules:
  //   MCtime must be inside the window [recoTime-deltaTimeMinus, recoTime+deltaTimePlus] to be associated to RecoHit
  //   MCtime must not already be in the window of the next hit, in which case it is associate to this one
  //
  //--------------------------------------------------------------------
  void CaloHitTruthMatch::makeTruthMatch(art::Event& event, CaloDigiMCCollection& caloDigiMCs,
                                         CaloDigiMCTruthAssn& caloDigiTruthMatch, CaloShowerMCTruthAssn& caloShowerTruthMatch)
  {
        
      art::ProductID digiMCProductID(event.getProductID<CaloDigiMCCollection>());
      const art::EDProductGetter* digiMCProductGetter = event.productGetter(digiMCProductID);

      const auto caloShowerSimHandle = event.getValidHandle(caloShowerSimToken_);
      const auto& caloShowerSims(*caloShowerSimHandle);
      const auto*  caloShowerSimBase  = &caloShowerSims.front();

      const auto caloCrystalHitHandle = event.getValidHandle(caloCrystalHitToken_);
      const auto& caloCrystalHits(*caloCrystalHitHandle);
      const auto* caloCrystalHitBase = &caloCrystalHits.front();
           

      std::map<int,std::vector<const CaloShowerSim*>> caloShowerSimsMap;
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
             std::vector<CaloEDepMC> edeps; 

             auto hitNextIt = std::next(hitIt);
             bool hitIsMatched(false);

             const CaloCrystalHit* hit = *hitIt;
             size_t idxHit             = (hit - caloCrystalHitBase);
             auto hitPtr               = art::Ptr<CaloCrystalHit>(caloCrystalHitHandle,idxHit);               

             if (diagLevel_ > 2) std::cout<<"[CaloHitTruthMatch]  inspect hit id/time/energy "<<hit->id()<<" / "<<hit->time()<<" / "<<hit->energyDep()<<std::endl;
                                      
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
                 const auto& sim    = showerSim->sim();
                 float edep         = showerSim->energyDep();
                 float edepG4       = showerSim->energyDepG4();
                 float time         = showerSim->time();
                 float pIn          = showerSim->momentumIn();
                 
                 //add shower to CaloEDepMC list
                 auto it = edeps.begin();
                 while (it != edeps.end()) {if (it->sim() == sim) break;  ++it;}

                 if (it!= edeps.end()) 
                    it->set(it->eDep()+edep, it->eDepG4()+edepG4,std::min(it->time(),time), std::max(pIn,it->momentumIn())); 
                 else                  
                    edeps.emplace_back(CaloEDepMC(sim,edep,edepG4,time,pIn));
                 
                 
                 // add shower to the detailed truthMatch
                 if (fillDetailedMC_) caloShowerTruthMatch.addSingle(hitPtr, showerSim->sim(), ShowerSimPtr);
                                 
                                                  
                 if (diagLevel_ > 1) diag(showerSim,hit);
                 if (diagLevel_ > 2) std::cout<<"[CaloHitTruthMatch]  matched shower id/time/energyDep()= "<<showerSim->crystalId()
                                              <<" / "<<showerSim->time()<<" / "<<showerSim->energyDep()<<std::endl;
                 ++showerIt;                
             }
             
             //sort CaloEDepMC by decreasing energy
             std::sort(edeps.begin(),edeps.end(),[](const auto& a, const auto& b){return a.eDep() > b.eDep();});

             caloDigiMCs.emplace_back(CaloDigiMC(std::move(edeps)));
             
             art::Ptr<CaloDigiMC> digiMCPtr = art::Ptr<CaloDigiMC>(digiMCProductID, caloDigiMCs.size(), digiMCProductGetter);             
             caloDigiTruthMatch.addSingle(hitPtr,digiMCPtr);
                          
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



