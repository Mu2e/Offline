//
// An EDProducer Module to match caloHits to MC info
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/MCDataProducts/inc/CaloEDepMC.hh"
#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include "Offline/MCDataProducts/inc/CaloShowerSim.hh"
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/Mu2eUtilities/inc/CaloPulseShape.hh"

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
             fhicl::Atom<art::InputTag>  caloHitCollection        { Name("caloHitCollection"),        Comment("Name of CaloHit collection") };
             fhicl::Atom<art::InputTag>  primaryParticle          { Name("primaryParticle"),          Comment("PrimaryParticle producer")};
             fhicl::Atom<std::string>    pulseFileName            { Name("pulseFileName"),            Comment("Calo pulse file name") };
             fhicl::Atom<std::string>    pulseHistName            { Name("pulseHistName"),            Comment("Calo pulse hist name") };
             fhicl::Atom<double>         digiSampling             { Name("digiSampling"),             Comment("Digitization time sampling") };
             fhicl::Atom<double>         minAmplitude             { Name("minAmplitude"),             Comment("Minimum amplitude of waveform to define hit length") };
             fhicl::Atom<double>         MeVToADC                 { Name("MeVToADC"),       Comment("MeV to ADC conversion factor") };
             fhicl::Atom<bool>           fillDetailedMC           { Name("fillDetailedMC"),           Comment("Fill SimParticle - SimShower Assn map")};
             fhicl::Atom<int>            diagLevel                { Name("diagLevel"),                Comment("Diag Level"),0 };
         };


        explicit CaloHitTruthMatch(const art::EDProducer::Table<Config>& config) :
           EDProducer{config},
           caloShowerSimToken_ {consumes<CaloShowerSimCollection> (config().caloShowerSimCollection())},
           caloHitToken_       {consumes<CaloHitCollection>(config().caloHitCollection())},
           ppToken_            {consumes<PrimaryParticle>(config().primaryParticle())},
           pulseFileName_      (config().pulseFileName()),
           pulseHistName_      (config().pulseHistName()),
           digiSampling_       (config().digiSampling()),
           minAmplitude_       (config().minAmplitude()),
           MeVToADC_           (config().MeVToADC()),
           fillDetailedMC_     (config().fillDetailedMC()),
           diagLevel_          (config().diagLevel())
        {
            produces<CaloHitMCCollection>();
            produces<CaloHitMCTruthAssn>();
            if (fillDetailedMC_) produces<CaloShowerMCTruthAssn>();
        }

        void beginJob() override;
        void produce(art::Event& e) override;
        void beginRun(art::Run& aRun) override;


      private:
         using SimParticlePtr = art::Ptr<SimParticle>;

         void makeTruthMatch (art::Event&, CaloHitMCCollection&, CaloHitMCTruthAssn&, CaloShowerMCTruthAssn&, const PrimaryParticle&);
         void fillEdeps      (const PrimaryParticle& primaryParticle, std::vector<CaloEDepMC>& edeps, const CaloShowerSim* showerSim);
         void diag           (const CaloShowerSim*, const CaloHit& );


         const art::ProductToken<CaloShowerSimCollection> caloShowerSimToken_;
         const art::ProductToken<CaloHitCollection>       caloHitToken_;
         const art::ProductToken<PrimaryParticle>         ppToken_;
         std::string                                      pulseFileName_;
         std::string                                      pulseHistName_;
         double                                           deltaTimeMinus_;
         double                                           digiSampling_;
         double                                           minAmplitude_;
         double                                           MeVToADC_;
         bool                                             fillDetailedMC_;
         std::vector<double>                              wf_;
         size_t                                           wfBinMax_;
         int                                              diagLevel_;

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

  //-----------------------------------------------------------------------------
  void CaloHitTruthMatch::beginRun(art::Run& aRun)
  {
      CaloPulseShape cps(pulseFileName_,pulseHistName_,digiSampling_);
      cps.buildShapes();

      wf_ = cps.digitizedPulse(0);
      for (auto& v : wf_) v *= MeVToADC_;
      wfBinMax_ = std::distance(wf_.begin(),std::max_element(wf_.begin(),wf_.end()));
  }


  //--------------------------------------------------------------------
  void CaloHitTruthMatch::produce(art::Event& event)
  {
      auto pph = event.getValidHandle<PrimaryParticle>(ppToken_);
      auto const& primaryParticles = *pph;

      std::unique_ptr<CaloHitMCCollection>   caloHitMCs(new CaloHitMCCollection);
      std::unique_ptr<CaloHitMCTruthAssn>    caloHitMCTruth(new CaloHitMCTruthAssn);
      std::unique_ptr<CaloShowerMCTruthAssn> caloShowerMCTruth(new CaloShowerMCTruthAssn);

      makeTruthMatch(event, *caloHitMCs, *caloHitMCTruth, *caloShowerMCTruth, primaryParticles);

      event.put(std::move(caloHitMCTruth));
      event.put(std::move(caloHitMCs));
      if (fillDetailedMC_) event.put(std::move(caloHitMCTruth));
  }



  // perform the association with the following rules:
  //   MCtime must be inside the window [recoTime-deltaTimeMinus, recoTime+deltaTimePlus] to be associated to RecoHit _unless_
  //   MCtime is already in the window of the next hit, in which case it is associate to this one.
  //--------------------------------------------------------------------
  void CaloHitTruthMatch::makeTruthMatch(art::Event& event, CaloHitMCCollection& caloHitMCs,
                                         CaloHitMCTruthAssn& CaloHitTruthMatch, CaloShowerMCTruthAssn& caloShowerTruthMatch,
                                         const PrimaryParticle& primaryParticle)
  {

      int nMatched(0);
      double totalEnergyMatched(0);

      // access collections of CaloHits and caloShowerHits, need art::ProductID for creating art::Ptr
      art::ProductID hitMCProductID(event.getProductID<CaloHitMCCollection>());
      const art::EDProductGetter* hitMCProductGetter = event.productGetter(hitMCProductID);
      const auto  caloHitHandle       = event.getValidHandle(caloHitToken_);
      const auto  caloShowerSimHandle = event.getValidHandle(caloShowerSimToken_);
      const auto& caloHits(*caloHitHandle);
      const auto& caloShowerSims(*caloShowerSimHandle);

      // sort the caloHits and caloShowerSim per crystal and then per time for each crystal to help with the matching algorithm.
      std::map<int, std::vector<const CaloHit*>>       caloHitMap;
      std::map<int, std::vector<const CaloShowerSim*>> caloShowerSimsMap;
      for (const auto& caloHit: caloHits) caloHitMap[caloHit.crystalID()].push_back(&caloHit);
      for (const auto& caloShowerSim: caloShowerSims) caloShowerSimsMap[caloShowerSim.crystalID()].push_back(&caloShowerSim);

      for (auto &kv : caloHitMap)        std::sort(kv.second.begin(),kv.second.end(), [](auto const a, auto const b){return a->time() < b->time();});
      for (auto &kv : caloShowerSimsMap) std::sort(kv.second.begin(),kv.second.end(), [](auto const a, auto const b){return a->time() < b->time();});


      // do the matching for a given calo hit
      for (size_t ihit=0; ihit < caloHits.size(); ++ihit)
      {
          const auto& hit        = caloHits[ihit];
          const auto sortedHits  = caloHitMap[hit.crystalID()];
          const auto sortedSims  = caloShowerSimsMap[hit.crystalID()];
          auto hitIt             = std::find(sortedHits.begin(),sortedHits.end(),&hit);
          auto hitNextIt         = std::next(hitIt);
          auto hitPtr            = art::Ptr<CaloHit>(caloHitHandle,ihit);
          auto showerIt          = sortedSims.begin();
          auto showerItEnd       = sortedSims.end();

          if (diagLevel_ > 2)
             for (const auto& shower : sortedSims) std::cout<<"[CaloHitTruthMatch] Sim shower  id/time/energy="<<shower->crystalID()
                                                          <<" / "<<shower->time()<<" / "<<shower->energyDep()<<std::endl;

          // Maximum time difference for a MChit to be associated to the reco hit given the reco hit amplitude and time of the next reco hit
          unsigned nbin(wfBinMax_);
          while (nbin < wf_.size()) {if (wf_[nbin]*hit.energyDep()<minAmplitude_) break; ++nbin;}
          double deltaTimePlus(nbin*digiSampling_);

          if (hitNextIt != sortedHits.end() && (*hitNextIt)->time()- (*hitIt)->time() - 2*digiSampling_ < deltaTimePlus )
          {
             deltaTimePlus = (*hitNextIt)->time()- (*hitIt)->time() - 2*digiSampling_;
          }
          if (diagLevel_ > 2) std::cout<<"[CaloHitTruthMatch]  inspect hit id/time/energy/length "<<hit.crystalID()
                                       <<" / "<<hit.time()<<" / "<<hit.energyDep()<<" "<<nbin*digiSampling_<<" "<<deltaTimePlus<<std::endl;


          //forward until we reach the recoHit time, then loop as long as the shower time is within the recoTime window
          bool hitIsMatched(false);
          std::vector<CaloEDepMC> edeps;
          while (showerIt != showerItEnd && ( (*showerIt)->time() < (*hitIt)->time() - 2*digiSampling_) ) ++showerIt;
          while (showerIt != showerItEnd && ( (*showerIt)->time() < (*hitIt)->time() + deltaTimePlus) )
          {
              hitIsMatched = true;
              const CaloShowerSim* showerSim = *showerIt;
              fillEdeps(primaryParticle, edeps, showerSim);

              if (fillDetailedMC_)
              {
                  size_t idxShower(0);
                  while (idxShower < caloShowerSims.size()) {if (&caloShowerSims[idxShower]==*showerIt) break; ++idxShower;}
                  auto  ShowerSimPtr = art::Ptr<CaloShowerSim>(caloShowerSimHandle,idxShower);
                  caloShowerTruthMatch.addSingle(hitPtr, showerSim->sim(), ShowerSimPtr);
              }

              if (diagLevel_ > 1) diag(showerSim,hit);
              if (diagLevel_ > 2) std::cout<<"[CaloHitTruthMatch]  matched shower id/time/energyDep()= "<<showerSim->crystalID()
                                           <<" / "<<showerSim->time()<<" / "<<showerSim->energyDep()<<std::endl;
              ++showerIt;
          }

          //sort CaloEDepMC by decreasing energy, create caloHit and keep track of hit -> MChit association
          std::sort(edeps.begin(),edeps.end(),[](const auto& a, const auto& b){return a.energyDep() > b.energyDep();});
          caloHitMCs.emplace_back(CaloHitMC(std::move(edeps)));

          art::Ptr<CaloHitMC> hitMCPtr = art::Ptr<CaloHitMC>(hitMCProductID, caloHitMCs.size()-1, hitMCProductGetter);
          CaloHitTruthMatch.addSingle(hitPtr,hitMCPtr);

          if (hitIsMatched) {totalEnergyMatched += (*hitIt)->energyDep();++nMatched;}
          ++hitIt;

          if (diagLevel_ > 2 && !hitIsMatched) std::cout<<"[CaloHitTruthMatch] hit not matched"<<std::endl;

      }

      if (diagLevel_ > 0) std::cout<<"[CaloHitTruthMatch]  total particles / energy matched = "<<nMatched<<" / "<<totalEnergyMatched<<std::endl;
  }


  //--------------------------------------------------------------------
  void CaloHitTruthMatch::fillEdeps(const PrimaryParticle& primaryParticle, std::vector<CaloEDepMC>& edeps, const CaloShowerSim* showerSim)
  {
      // check if there is already a caloEdep object with same caloShowerSim's SimParticle
      auto it = edeps.begin();
      while (it != edeps.end()) {if (it->sim() == showerSim->sim()) break;  ++it;}

      // if found the add caloShowerSim to matching caloEdep object, otherwise create new caloEdep
      if (it!= edeps.end())
      {
          it->addEDep(showerSim->energyDep());
          it->addEDepG4(showerSim->energyDepG4());
          it->addTime(showerSim->time());
          it->addMom(showerSim->momentumIn());
      }
      else
      {
          MCRelationship mcrel;
          for (const auto& spp : primaryParticle.primarySimParticles())
          {
              MCRelationship mcr(spp,showerSim->sim());
              if (mcr > mcrel) mcrel = mcr;
          }
          edeps.emplace_back(CaloEDepMC(showerSim->sim(),showerSim->energyDep(),showerSim->energyDepG4(),
                                        showerSim->time(),showerSim->momentumIn(),mcrel));
      }
  }



  //--------------------------------------------------------------------
  void CaloHitTruthMatch::diag(const CaloShowerSim* shower, const CaloHit& hit)
  {
      hTime_->Fill(shower->time()-hit.time());
      hEnerTime_->Fill(shower->energyDep(),shower->time()-hit.time());
      hTime2d_->Fill(shower->time(),hit.time());
      hEner2d_->Fill(shower->energyDep(),hit.energyDep());
      hdEdT_->Fill(hit.energyDep()-shower->energyDep(),shower->time()-hit.time());

      if (shower->energyDep() > 5) hTime2_->Fill(shower->time()-hit.time());

      double deltaE = std::abs(shower->energyDep()-hit.energyDep());
      if (deltaE > 5 && shower->energyDep() > 5)
      if (!hit.recoCaloDigis().empty()) hChi2_->Fill(hit.recoCaloDigis().at(0)->chi2()/hit.recoCaloDigis().at(0)->ndf());
  }


}

using mu2e::CaloHitTruthMatch;
DEFINE_ART_MODULE(CaloHitTruthMatch)



