//
// An EDProducer Module to match caloHits to MC info
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/types/Atom.h"

#include "Offline/MCDataProducts/inc/CaloEDepMC.hh"
#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include "Offline/MCDataProducts/inc/CaloShowerSim.hh"
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/MCRelationship.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/Mu2eUtilities/inc/CaloPulseUtil.hh"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>


namespace mu2e {

  class CaloHitTruthMatch : public art::EDProducer
  {
      public:
         struct Config
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             using CPG     = CaloPulseUtil::Config;
             fhicl::Table<CPG>          pulseCache              { Name("pulseCache"),              Comment("Pulse cache maker config") };
             fhicl::Atom<art::InputTag> caloShowerSimCollection { Name("caloShowerSimCollection"), Comment("Name of caloShowerSim Collection") };
             fhicl::Atom<art::InputTag> caloHitCollection       { Name("caloHitCollection"),       Comment("Name of CaloHit collection") };
             fhicl::Atom<art::InputTag> primaryParticle         { Name("primaryParticle"),         Comment("PrimaryParticle producer")};
             fhicl::Atom<double>        digiSampling            { Name("digiSampling"),            Comment("Digitizer sampling time (ns) ") };
             fhicl::Atom<double>        deltaTimeMinus          { Name("deltaTimeMinus"),          Comment("Max time (ns) a MC hit may precede the reco hit to be matched"), 100.0 };
             fhicl::Atom<double>        minAmplitude            { Name("minAmplitude"),            Comment("Minimum amplitude of waveform to define hit length") };
             fhicl::Atom<bool>          fillDetailedMC          { Name("fillDetailedMC"),          Comment("Fill SimParticle - SimShower Assn map")};
             fhicl::Atom<int>           diagLevel               { Name("diagLevel"),               Comment("Diag Level"),0 };
         };


        explicit CaloHitTruthMatch(const art::EDProducer::Table<Config>& config) :
           EDProducer{config},
           caloShowerSimToken_ {consumes<CaloShowerSimCollection>(config().caloShowerSimCollection())},
           caloHitToken_       {consumes<CaloHitCollection>(config().caloHitCollection())},
           ppToken_            {consumes<PrimaryParticle>(config().primaryParticle())},
           pulseCache_         (config().pulseCache()),
           digiSampling_       (config().digiSampling()),
           deltaTimeMinus_     (config().deltaTimeMinus()),
           minAmplitude_       (config().minAmplitude()),
           fillDetailedMC_     (config().fillDetailedMC()),
           diagLevel_          (config().diagLevel())
        {
            produces<CaloHitMCCollection>();
            produces<CaloHitMCTruthAssn>();
            if (fillDetailedMC_) produces<CaloShowerMCTruthAssn>();
        }

        void beginRun(art::Run&) override;
        void produce (art::Event& e) override;


      private:
         using SimParticlePtr = art::Ptr<SimParticle>;

         void makeTruthMatch(art::Event&, CaloHitMCCollection&, CaloHitMCTruthAssn&, CaloShowerMCTruthAssn&, const PrimaryParticle&);
         void fillEdeps     (const PrimaryParticle& primaryParticle, std::vector<CaloEDepMC>& edeps, const CaloShowerSim* showerSim);

         const art::ProductToken<CaloShowerSimCollection> caloShowerSimToken_;
         const art::ProductToken<CaloHitCollection>       caloHitToken_;
         const art::ProductToken<PrimaryParticle>         ppToken_;
         CaloPulseUtil::Config                            pulseCache_;
         double                                           digiSampling_;
         double                                           deltaTimeMinus_;
         double                                           minAmplitude_;
         bool                                             fillDetailedMC_;
         std::vector<double>                              wf_;
         std::size_t                                      wfBinMax_{0};
         int                                              diagLevel_;
   };




  //-----------------------------------------------------------------------------
  void CaloHitTruthMatch::beginRun(art::Run&)
  {
      CaloPulseUtil cps(pulseCache_);
      cps.buildCache();

      wf_       = cps.digitizedPulse(0);
      wfBinMax_ = std::distance(wf_.begin(),std::max_element(wf_.begin(),wf_.end()));
  }


  //--------------------------------------------------------------------
  void CaloHitTruthMatch::produce(art::Event& event)
  {
      const auto& primaryParticles = *event.getValidHandle<PrimaryParticle>(ppToken_);

      auto caloHitMCs        = std::make_unique<CaloHitMCCollection>();
      auto caloHitMCTruth    = std::make_unique<CaloHitMCTruthAssn>();
      auto caloShowerMCTruth = std::make_unique<CaloShowerMCTruthAssn>();

      makeTruthMatch(event, *caloHitMCs, *caloHitMCTruth, *caloShowerMCTruth, primaryParticles);

      event.put(std::move(caloHitMCs));
      event.put(std::move(caloHitMCTruth));
      if (fillDetailedMC_) event.put(std::move(caloShowerMCTruth));
  }



  // Association rules:
  //   MCtime must be inside [recoTime-deltaTimeMinus, recoTime+deltaTimePlus] to be associated to a RecoHit,
  //   _unless_ MCtime is already inside the window of the next hit, in which case it goes to that one.
  //--------------------------------------------------------------------
  void CaloHitTruthMatch::makeTruthMatch(art::Event& event, CaloHitMCCollection& caloHitMCs,
                                         CaloHitMCTruthAssn& caloHitTruthMatch, CaloShowerMCTruthAssn& caloShowerTruthMatch,
                                         const PrimaryParticle& primaryParticle)
  {
      int    nMatched(0);
      double totalEnergyMatched(0);

      const art::ProductID hitMCProductID(event.getProductID<CaloHitMCCollection>());
      const art::EDProductGetter* hitMCProductGetter = event.productGetter(hitMCProductID);

      const auto  caloHitHandle       = event.getValidHandle(caloHitToken_);
      const auto  caloShowerSimHandle = event.getValidHandle(caloShowerSimToken_);
      const auto& caloHits(*caloHitHandle);
      const auto& caloShowerSims(*caloShowerSimHandle);

      // sort caloHits and caloShowerSims per crystal, then per time, to help the matching
      std::unordered_map<int, std::vector<const CaloHit*>>       caloHitMap;
      std::unordered_map<int, std::vector<const CaloShowerSim*>> caloShowerSimsMap;
      for (const auto& caloHit       : caloHits)       caloHitMap[caloHit.crystalID()].push_back(&caloHit);
      for (const auto& caloShowerSim : caloShowerSims) caloShowerSimsMap[caloShowerSim.crystalID()].push_back(&caloShowerSim);

      const auto byTime = [](const auto* a, const auto* b){return a->time() < b->time();};
      for (auto& [id, v] : caloHitMap)        std::sort(v.begin(),v.end(),byTime);
      for (auto& [id, v] : caloShowerSimsMap) std::sort(v.begin(),v.end(),byTime);

      // O(1) shower-pointer -> index lookup for the detailed-MC association (was an O(n) scan per match)
      std::unordered_map<const CaloShowerSim*,std::size_t> showerIndex;
      if (fillDetailedMC_) {
          showerIndex.reserve(caloShowerSims.size());
          for (std::size_t i=0; i<caloShowerSims.size(); ++i) showerIndex.emplace(&caloShowerSims[i], i);
      }

      // matching per calo hit
      caloHitMCs.reserve(caloHits.size());
      for (std::size_t ihit=0; ihit < caloHits.size(); ++ihit) {
          const CaloHit& hit          = caloHits[ihit];
          const auto&    sortedHits   = caloHitMap.at(hit.crystalID());
          const auto&    sortedSims   = caloShowerSimsMap.at(hit.crystalID());
          const auto     hitIt        = std::find(sortedHits.begin(),sortedHits.end(),&hit);
          const auto     hitNextIt    = std::next(hitIt);
          const auto     hitPtr       = art::Ptr<CaloHit>(caloHitHandle,ihit);

          if (diagLevel_ > 2)
             for (const auto* shower : sortedSims)
                std::cout<<"[CaloHitTruthMatch] Sim shower id/time/energy="<<shower->crystalID()
                         <<" / "<<shower->time()<<" / "<<shower->energyDep()<<std::endl;

          // Maximum time difference for an MC hit to be associated, given the reco hit amplitude and the next reco hit time
          std::size_t nbin(wfBinMax_);
          while (nbin < wf_.size()) {if (wf_[nbin]*hit.energyDep() < minAmplitude_) break; ++nbin;}
          double deltaTimePlus(nbin*digiSampling_);

          if (hitNextIt != sortedHits.end() && (*hitNextIt)->time() - hit.time() - 2*digiSampling_ < deltaTimePlus)
             deltaTimePlus = (*hitNextIt)->time() - hit.time() - 2*digiSampling_;

          if (diagLevel_ > 2) std::cout<<"[CaloHitTruthMatch] inspect hit id/time/energy/length "<<hit.crystalID()
                                       <<" / "<<hit.time()<<" / "<<hit.energyDep()<<" "<<nbin*digiSampling_<<" "<<deltaTimePlus<<std::endl;

          // forward to the reco hit time, then loop while the shower time is within the reco window
          bool hitIsMatched(false);
          std::vector<CaloEDepMC> edeps;
          auto showerIt = sortedSims.begin();
          while (showerIt != sortedSims.end() && (*showerIt)->time() < hit.time() - deltaTimeMinus_) ++showerIt;
          while (showerIt != sortedSims.end() && (*showerIt)->time() < hit.time() + deltaTimePlus) {
              hitIsMatched = true;
              const CaloShowerSim* showerSim = *showerIt;
              fillEdeps(primaryParticle, edeps, showerSim);

              if (fillDetailedMC_) {
                  const auto showerSimPtr = art::Ptr<CaloShowerSim>(caloShowerSimHandle, showerIndex.at(showerSim));
                  caloShowerTruthMatch.addSingle(hitPtr, showerSim->sim(), showerSimPtr);
              }

              if (diagLevel_ > 2) std::cout<<"[CaloHitTruthMatch] matched shower id/time/energyDep= "<<showerSim->crystalID()
                                           <<" / "<<showerSim->time()<<" / "<<showerSim->energyDep()<<std::endl;
              ++showerIt;
          }

          //sort CaloEDepMC by decreasing energy, create caloHit and keep track of hit -> MChit association
          std::sort(edeps.begin(),edeps.end(),[](const auto& a, const auto& b){return a.energyDep() > b.energyDep();});
          caloHitMCs.emplace_back(std::move(edeps),hit.crystalID());

          const auto hitMCPtr = art::Ptr<CaloHitMC>(hitMCProductID, caloHitMCs.size()-1, hitMCProductGetter);
          caloHitTruthMatch.addSingle(hitPtr,hitMCPtr);

          if (hitIsMatched) {totalEnergyMatched += hit.energyDep(); ++nMatched;}
          else if (diagLevel_ > 2) std::cout<<"[CaloHitTruthMatch] hit not matched"<<std::endl;
      }

      if (diagLevel_ > 0) std::cout<<"[CaloHitTruthMatch] total particles / energy matched = "<<nMatched<<" / "<<totalEnergyMatched<<std::endl;
  }


  //--------------------------------------------------------------------
  void CaloHitTruthMatch::fillEdeps(const PrimaryParticle& primaryParticle, std::vector<CaloEDepMC>& edeps, const CaloShowerSim* showerSim)
  {
      // if there is already a CaloEDepMC for this SimParticle, merge; otherwise create a new one
      auto it = std::find_if(edeps.begin(), edeps.end(),
                             [&](const CaloEDepMC& e){return e.sim() == showerSim->sim();});

      if (it != edeps.end()) {
          it->addEDep  (showerSim->energyDep());
          it->addEDepG4(showerSim->energyDepG4());
          it->addTime  (showerSim->time());
          it->addMom   (showerSim->momentumIn());
      }
      else{
          MCRelationship mcrel;
          for (const auto& spp : primaryParticle.primarySimParticles()) {
              MCRelationship mcr(spp,showerSim->sim());
              if (mcr > mcrel) mcrel = mcr;
          }
          edeps.emplace_back(showerSim->sim(),showerSim->energyDep(),showerSim->energyDepG4(),
                             showerSim->time(),showerSim->momentumIn(),mcrel);
      }
  }

}

DEFINE_ART_MODULE(mu2e::CaloHitTruthMatch)
