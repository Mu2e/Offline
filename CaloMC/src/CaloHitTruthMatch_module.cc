//
// An EDProducer Module to match caloCrystalHits to MC info
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "fhiclcpp/types/Atom.h"

#include "MCDataProducts/inc/CaloEDepMC.hh"
#include "MCDataProducts/inc/CaloDigiMC.hh"
#include "MCDataProducts/inc/CaloShowerStep.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "MCDataProducts/inc/PrimaryParticle.hh"
#include "MCDataProducts/inc/SimParticle.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"

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
             fhicl::Atom<art::InputTag>  primaryParticle          { Name("primaryParticle"),	      Comment("PrimaryParticle producer")};
             fhicl::Atom<double>         deltaTimeMinus           { Name("deltaTimeMinus"),           Comment("Maximum time before hit to include MC hit") };
             fhicl::Atom<double>         deltaTimePlus            { Name("deltaTimePlus"),            Comment("Maximum time after hit to include MC hit") };
             fhicl::Atom<bool>           fillDetailedMC           { Name("fillDetailedMC"),           Comment("Fill SimParticle - SimShower Assn map")};
             fhicl::Atom<int>            diagLevel                { Name("diagLevel"),                Comment("Diag Level"),0 };
         };


        explicit CaloHitTruthMatch(const art::EDProducer::Table<Config>& config) :
           EDProducer{config},
           caloShowerSimToken_ {consumes<CaloShowerSimCollection> (config().caloShowerSimCollection())},
           caloCrystalHitToken_{consumes<CaloCrystalHitCollection>(config().caloCrystalHitCollection())},
           ppToken_            {consumes<PrimaryParticle>(config().primaryParticle())},
           deltaTimeMinus_     (config().deltaTimeMinus()),
           deltaTimePlus_      (config().deltaTimePlus()),
           fillDetailedMC_     (config().fillDetailedMC()),
           diagLevel_          (config().diagLevel())
        {
            produces<CaloDigiMCCollection>();    
            produces<CaloDigiMCTruthAssn>();    
            if (fillDetailedMC_) produces<CaloShowerMCTruthAssn>();    
        }

        void beginJob() override;
        void produce(art::Event& e) override;


      private:
         using SimParticlePtr = art::Ptr<SimParticle>;

         void makeTruthMatch (art::Event&, CaloDigiMCCollection&, CaloDigiMCTruthAssn&, CaloShowerMCTruthAssn&);
         void diag           (const CaloShowerSim*, const CaloCrystalHit* );


         const art::ProductToken<CaloShowerSimCollection>  caloShowerSimToken_;
         const art::ProductToken<CaloCrystalHitCollection> caloCrystalHitToken_;
         const art::ProductToken<PrimaryParticle>          ppToken_;
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

      const auto caloCrystalHitHandle = event.getValidHandle(caloCrystalHitToken_);
      const auto& caloCrystalHits(*caloCrystalHitHandle);
      const auto caloShowerSimHandle = event.getValidHandle(caloShowerSimToken_);
      const auto& caloShowerSims(*caloShowerSimHandle);
      //const auto primaryParticleHandle = event.getValidHandle<PrimaryParticle>(ppToken_);
      //const auto& primaryParticle = *primaryParticleHandle;
     
      
      std::map<int, std::vector<const CaloCrystalHit*>> caloHitMap;
      std::map<int,std::vector<const CaloShowerSim*>>   caloShowerSimsMap;
      for (auto const& CaloCrystalHit: caloCrystalHits) caloHitMap[CaloCrystalHit.id()].push_back(&CaloCrystalHit);
      for (auto const& caloShowerSim: caloShowerSims)   caloShowerSimsMap[caloShowerSim.crystalId()].push_back(&caloShowerSim);
      
      
      int nMatched(0);
      double totalEnergyMatched(0); 
      for (auto &kv : caloHitMap)
      {
          int crystalId = kv.first;

          std::vector<const CaloCrystalHit*> &caloHits = kv.second;          
          std::sort(caloHits.begin(),caloHits.end(), [](auto const a, auto const b){return a->time() < b->time();});
          
          std::vector<const CaloShowerSim*>& caloSims = caloShowerSimsMap[crystalId];
          std::sort(caloSims.begin(),caloSims.end(), [](auto const a, auto const b){return a->time() < b->time();});
          
          if (diagLevel_ > 2) 
             for (const auto& shower : caloSims) 
                 std::cout<<"[CaloHitTruthMatch] Sim shower  id/energy/time="<<shower->crystalId()<<" / "<<shower->energyDep()<<" / "<<shower->time()<<std::endl;
          

          auto showerIt    = caloSims.begin();
          auto showerItEnd = caloSims.end();
          auto hitIt       = caloHits.begin();
          auto hitItEnd    = caloHits.end();

          while (hitIt != hitItEnd)
          {            
             std::vector<CaloEDepMC> edeps; 

             auto hitNextIt = std::next(hitIt);
             bool hitIsMatched(false);
             
             const CaloCrystalHit* hit = *hitIt;
             size_t idxHit(0); 
             while (idxHit < caloCrystalHits.size()) {if (&caloCrystalHits[idxHit]==hit) break; ++idxHit;}
             auto hitPtr = art::Ptr<CaloCrystalHit>(caloCrystalHitHandle,idxHit);               

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
                 
                 size_t idxShower(0); 
                 while (idxShower < caloShowerSims.size()) {if (&caloShowerSims[idxShower]==*showerIt) break; ++idxShower;}
                 auto  ShowerSimPtr = art::Ptr<CaloShowerSim>(caloShowerSimHandle,idxShower);
                 
                 const CaloShowerSim* showerSim(*showerIt);                 
                 const auto& sim    = showerSim->sim();
                 float edep         = showerSim->energyDep();
                 float edepG4       = showerSim->energyDepG4();
                 float time         = showerSim->time();
                 float pIn          = showerSim->momentumIn();
                 
                 //add shower to CaloEDepMC list
                 auto it = edeps.begin();
                 while (it != edeps.end()) {if (it->sim() == sim) break;  ++it;}

                 if (it!= edeps.end())
                 { 
                    it->addEDep(edep);
                    it->addEDepG4(edepG4);
                    it->addTime(time);
                    it->addMom(pIn);
                 }
                 else 
                 {                                     
                    MCRelationship mcrel;
                    //for (const auto& spp : primaryParticle.primarySimParticles())
                    //{
                    //    MCRelationship mcr(spp,sim);
                    //    if (mcr > mcrel) mcrel = mcr;
                    //}
                    edeps.emplace_back(CaloEDepMC(sim,edep,edepG4,time,pIn,mcrel));
                 }
                 
                 // add shower to the detailed truthMatch
                 if (fillDetailedMC_) caloShowerTruthMatch.addSingle(hitPtr, showerSim->sim(), ShowerSimPtr);
                                 
                                                  
                 if (diagLevel_ > 1) diag(showerSim,hit);
                 if (diagLevel_ > 2) std::cout<<"[CaloHitTruthMatch]  matched shower id/time/energyDep()= "<<showerSim->crystalId()
                                              <<" / "<<showerSim->time()<<" / "<<showerSim->energyDep()<<std::endl;
                 ++showerIt;                
             }
             
             //sort CaloEDepMC by decreasing energy
             std::sort(edeps.begin(),edeps.end(),[](const auto& a, const auto& b){return a.energyDep() > b.energyDep();});
             caloDigiMCs.emplace_back(CaloDigiMC(std::move(edeps)));

             art::Ptr<CaloDigiMC> digiMCPtr = art::Ptr<CaloDigiMC>(digiMCProductID, caloDigiMCs.size()-1, digiMCProductGetter);             
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



