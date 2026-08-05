//
// Transform the energy deposited in the scintillator into photo-electrons (PE) seen by the photosensor.
// Includes corrections from Birks law, longitudinal response uniformity and photo-statistics fluctuations.
// The PE are generated individually and corrected for transit time.
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/CaloMC/inc/CaloPhotonPropagation.hh"
#include "Offline/CaloConditions/inc/CalSimParams.hh"
#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/DataProducts/inc/CrystalId.hh"
#include "Offline/ProditionsService/inc/ProditionsHandle.hh"
#include "Offline/DataProducts/inc/EventWindowMarker.hh"
#include "Offline/DAQConditions/inc/EventTiming.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/PhysicsParams.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CaloShowerRO.hh"
#include "Offline/MCDataProducts/inc/CaloShowerSim.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/SeedService/inc/SeedService.hh"
#include "CLHEP/Random/RandPoissonQ.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <set>
#include <utility>
#include <vector>


namespace {

  struct StepEntry
  {
      StepEntry(const art::Ptr<mu2e::CaloShowerStep>& step, float edepCorr, float timeCorr) :
        step_(step),edepCorr_(edepCorr),timeCorr_(timeCorr)
      {}
      art::Ptr<mu2e::CaloShowerStep> step_;
      float edepCorr_,timeCorr_;
  };

  struct SimParticleSummary
  {
      SimParticleSummary(const art::Ptr<mu2e::CaloShowerStep>& step, float edepCorr, float timeCorr) :
        steps_{step},edepCorr_(edepCorr),timeCorr_(timeCorr)
      {}
      void add(const art::Ptr<mu2e::CaloShowerStep>& step, float edepCorr, float timeCorr)
      {
         steps_.push_back(step);
         edepCorr_ += edepCorr;
         timeCorr_  = std::min(timeCorr,timeCorr_);
      }
      std::vector<art::Ptr<mu2e::CaloShowerStep>> steps_;
      float edepCorr_,timeCorr_;
  };

  struct DiagSummary
  {
      int   totSteps{0}, totNPE{0};
      float totEdep{0.f}, totEdepCorr{0.f}, totEdepNPE{0.f};
  };
}



namespace mu2e {

  class CaloShowerROMaker : public art::EDProducer
  {
      public:
         struct Config
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             fhicl::Sequence<art::InputTag>  caloShowerStepCollection { Name("caloShowerStepCollection"), Comment("Compressed shower inputs for calo crystals") };
             fhicl::Atom<art::InputTag>      ewMarkerTag              { Name("eventWindowMarker"),        Comment("EventWindowMarker producer") };
             fhicl::Atom<art::InputTag>      pbtmcTag                 { Name("protonBunchTimeMC"),        Comment("ProtonBunchTimeMC producer") };
             fhicl::Atom<std::string>        propagationFileName      { Name("propagationFileName"),      Comment("Photon propagation file name")};
             fhicl::Atom<std::string>        propagationHistName      { Name("propagationHistName"),      Comment("Photon propagation hist name")};
             fhicl::Atom<float>              digitizationStart        { Name("digitizationStart"),        Comment("Minimum time of hit to be digitized") };
             fhicl::Atom<float>              digitizationEnd          { Name("digitizationEnd"),          Comment("Maximum time of hit to be digitized") };
             fhicl::Atom<float>              digitizationBuffer       { Name("digitizationBuffer"),       Comment("Digi time buffer for photon propagation inside crystal") };
             fhicl::Atom<bool>               LRUCorrection            { Name("LRUCorrection"),            Comment("Include LRU corrections") };
             fhicl::Atom<bool>               BirksCorrection          { Name("BirksCorrection"),          Comment("Include Birks corrections") };
             fhicl::Atom<bool>               PEStatCorrection         { Name("PEStatCorrection"),         Comment("Include PE Poisson fluctuations") };
             fhicl::Atom<bool>               addTravelTime            { Name("addTravelTime"),            Comment("Include light propagation time") };
             fhicl::Atom<int>                diagLevel                { Name("diagLevel"),                Comment("Diag Level"),0 };
         };

         explicit CaloShowerROMaker(const art::EDProducer::Table<Config>& config) :
            EDProducer{config},
            ewMarkerToken_       {consumes<EventWindowMarker>(config().ewMarkerTag())},
            pbtmcToken_          {consumes<ProtonBunchTimeMC>(config().pbtmcTag())},
            digitizationStart_   (config().digitizationStart()),
            digitizationEnd_     (config().digitizationEnd()),
            digitizationBuffer_  (config().digitizationBuffer()),
            LRUCorrection_       (config().LRUCorrection()),
            BirksCorrection_     (config().BirksCorrection()),
            PEStatCorrection_    (config().PEStatCorrection()),
            addTravelTime_       (config().addTravelTime()),
            diagLevel_           (config().diagLevel()),
            engine_              (createEngine(art::ServiceHandle<SeedService>()->getSeed())),
            randPoisson_         (engine_),
            photonProp_          (config().propagationFileName(),config().propagationHistName(),engine_)
         {
             for (const auto& tag : config().caloShowerStepCollection())
                crystalShowerTokens_.push_back(consumes<CaloShowerStepCollection>(tag));

             produces<CaloShowerROCollection>();
             produces<CaloShowerSimCollection>();
         }

         void beginRun(art::Run& aRun) override;
         void produce(art::Event& e)  override;


      private:
         using StepHandles = std::vector<art::ValidHandle<CaloShowerStepCollection>>;

         void  makeReadoutHits (const StepHandles&, CaloShowerROCollection&, CaloShowerSimCollection&, const CalSimParams&,
                                const EventWindowMarker&, const ProtonBunchTimeMC&, float timeFromProtonsToDRMarker);
         float LRUCorrection   (float normalizedPosZ, float edepInit, float lru) const;
         void  dumpCaloShowerSim(const CaloShowerSimCollection& caloShowerSims) const;

         ProditionsHandle<CalSimParams>                          calCrystalConds_;
         std::vector<art::ProductToken<CaloShowerStepCollection>> crystalShowerTokens_;
         art::ProductToken<EventWindowMarker> ewMarkerToken_;
         art::ProductToken<ProtonBunchTimeMC> pbtmcToken_;
         float                   digitizationStart_;
         float                   digitizationEnd_;
         float                   digitizationBuffer_;
         bool                    LRUCorrection_;
         bool                    BirksCorrection_;
         bool                    PEStatCorrection_;
         bool                    addTravelTime_;
         int                     diagLevel_;
         CLHEP::HepRandomEngine& engine_;
         CLHEP::RandPoissonQ     randPoisson_;
         CaloPhotonPropagation   photonProp_;
  };



  //-----------------------------------------------
  void CaloShowerROMaker::beginRun(art::Run&)
  {
      photonProp_.buildTable();
  }


  //---------------------------------------------------------------
  void CaloShowerROMaker::produce(art::Event& event)
  {
      if (diagLevel_ > 0) std::cout << "[CaloShowerROMaker::produce] begin" << std::endl;

      const EventWindowMarker& ewMarker = *event.getValidHandle(ewMarkerToken_);
      const ProtonBunchTimeMC& pbtmc    = *event.getValidHandle(pbtmcToken_);

      ProditionsHandle<EventTiming> eventTimingHandle;
      const EventTiming& eventTiming = eventTimingHandle.get(event.id());
      const float timeFromProtonsToDRMarker = eventTiming.timeFromProtonsToDRMarker();

      const auto& calCrystalConds = calCrystalConds_.get(event.id());

      auto caloShowerROs  = std::make_unique<CaloShowerROCollection>();
      auto caloShowerSims = std::make_unique<CaloShowerSimCollection>();

      StepHandles crystalShowerHandles;
      crystalShowerHandles.reserve(crystalShowerTokens_.size());
      std::transform(crystalShowerTokens_.begin(), crystalShowerTokens_.end(), std::back_inserter(crystalShowerHandles),
                     [&event](const auto& token){return event.getValidHandle(token);});

      makeReadoutHits(crystalShowerHandles, *caloShowerROs, *caloShowerSims, calCrystalConds, ewMarker, pbtmc, timeFromProtonsToDRMarker);

      event.put(std::move(caloShowerROs));
      event.put(std::move(caloShowerSims));

      if (diagLevel_ > 0) std::cout << "[CaloShowerROMaker::produce] end" << std::endl;
  }


  //-----------------------------------------------------------------------------------------------------
  void CaloShowerROMaker::makeReadoutHits(const StepHandles& crystalShowerHandles, CaloShowerROCollection& caloShowerROs,
                                          CaloShowerSimCollection& caloShowerSims, const CalSimParams& calCrystalConds,
                                          const EventWindowMarker& ewMarker, const ProtonBunchTimeMC& pbtmc,
                                          float timeFromProtonsToDRMarker)
  {
      GlobalConstantsHandle<ParticleDataList> pdt;
      const float mbtime = GlobalConstantsHandle<PhysicsParams>()->getNominalDRPeriod();

      const Calorimeter& cal    = *(GeomHandle<Calorimeter>());
      const float crystalLength = cal.G4Info().get<double>("crystalZLength");

      std::map<int,std::vector<StepEntry>> simEntriesMap;
      DiagSummary diagSum;

      // Digitization start / end from accelerator DR marker with PB jitter
      const float correctedDigitizeStart = digitizationStart_ - pbtmc.pbtime_ - timeFromProtonsToDRMarker - digitizationBuffer_;
      const float correctedDigitizeEnd   = digitizationEnd_   - pbtmc.pbtime_ - timeFromProtonsToDRMarker;

      //-----------------------------------------------------------------------
      // store corrected energy deposits for each readout
      for (const auto& showerHandle : crystalShowerHandles)
      {
          const CaloShowerStepCollection& caloShowerSteps(*showerHandle);
          for (auto istep = caloShowerSteps.begin(); istep != caloShowerSteps.end(); ++istep)
          {
              const CaloShowerStep& step = *istep;

              // see doc-db for calo folding description. Note pbtmc.pbtime_ is NEGATIVE!
              const double mbLength   = (ewMarker.spillType() == EventWindowMarker::SpillType::onspill) ? mbtime : ewMarker.eventLength();
              const double maxHitTime = std::max<double>(mbtime,correctedDigitizeEnd);

              double hitTime = std::fmod(step.time()+pbtmc.pbtime_, mbLength) - pbtmc.pbtime_;
              if (hitTime < maxHitTime-mbLength) hitTime += mbLength;
              if (hitTime < correctedDigitizeStart) continue;

              const std::size_t idx = std::distance(caloShowerSteps.begin(), istep);
              const auto stepPtr = art::Ptr<CaloShowerStep>(showerHandle,idx);

              const auto  crystalID  = CrystalId(step.volumeG4ID());
              const int   SiPMIDBase = crystalID.SiPMId(CaloConst::SiPM0);
              const float posZ       = step.position().z();
              const float lru        = calCrystalConds.LRU(crystalID);
              const auto  pePerMeVs  = calCrystalConds.pePerMeVs(crystalID);

              float edep_corr(step.energyDepG4());
              if (BirksCorrection_) edep_corr = step.energyDepBirks();
              if (LRUCorrection_)   edep_corr = LRUCorrection(posZ/crystalLength, edep_corr, lru);

              // Generate individual PEs and their arrival times
              for (int i=0; i<CaloConst::_nSiPMPerCrystal; ++i)
              {
                  const int   SiPMID = SiPMIDBase + i;
                  const float mean   = std::max(0.0f, edep_corr*pePerMeVs[i]);   // guard: negative mean is UB in RandPoissonQ
                  const int   NPE    = PEStatCorrection_ ? randPoisson_.fire(mean) : int(std::lround(mean));
                  if (NPE <= 0) continue;

                  std::vector<float> PETime(NPE,hitTime);
                  if (addTravelTime_)
                     for (auto& t : PETime) t += photonProp_.propTimeSimu(crystalLength-posZ);

                  if (diagLevel_ > 2)
                  {
                      std::cout<<"[CaloShowerROMaker] SiPMID:"<<SiPMID<<"  energy / NPE = "<<edep_corr<<"  /  "<<NPE<<"\nTime hit "<<std::endl;
                      for (float t : PETime) std::cout<<t<<" ";
                      std::cout<<std::endl;
                  }

                  diagSum.totNPE     += NPE;
                  diagSum.totEdepNPE += double(NPE)/pePerMeVs[i]/2.0;   // average between the two RO

                  caloShowerROs.emplace_back(SiPMID, stepPtr, std::move(PETime));
              }
              diagSum.totEdep     += step.energyDepG4();
              diagSum.totEdepCorr += edep_corr;
              diagSum.totSteps    += step.nCompress();

              simEntriesMap[crystalID.id()].emplace_back(stepPtr,edep_corr,hitTime);
          }
      }

      // sort once, after all input handles are processed
      std::sort(caloShowerROs.begin(),caloShowerROs.end(),
                [](const auto& a, const auto& b){return a.SiPMID() < b.SiPMID();});

      //--------------------------------------------------
      // Produce the final MC truth info collecting energy deposits for each SimParticle in each crystal
      for (const auto& [crId, newSteps] : simEntriesMap)
      {
          std::map<art::Ptr<SimParticle>,SimParticleSummary> summaryMap;
          for (const auto& newStep : newSteps)
          {
              const art::Ptr<SimParticle>& sim = newStep.step_->simParticle();
              auto mfind = summaryMap.find(sim);
              if (mfind==summaryMap.end())
                 summaryMap.emplace(sim, SimParticleSummary(newStep.step_,newStep.edepCorr_,newStep.timeCorr_));
              else
                 mfind->second.add(newStep.step_,newStep.edepCorr_,newStep.timeCorr_);
          }

          for (auto& [sim, summ] : summaryMap)
             caloShowerSims.emplace_back(summ.steps_, summ.edepCorr_, summ.timeCorr_);
      }

      //--------------------------------------------------
      // Diag
      if (diagLevel_ > 2) dumpCaloShowerSim(caloShowerSims);

      if (diagLevel_ > 1)
      {
         std::set<int> crIds;
         for (const auto& css : caloShowerSims) crIds.insert(css.crystalID());

         for (int crId : crIds)
         {
            std::map<art::Ptr<SimParticle>, double> simMap;
            for (const auto& css : caloShowerSims) if (css.crystalID()==crId) simMap[css.sim()] += css.energyDep();
            for (const auto& [sim, energy] : simMap) std::cout<<"CrId: "<<crId<<"  Sim id: "<<sim.id()<<"   energy="<<energy<<std::endl;
         }
      }

      if (diagLevel_ > 0) std::cout<<"[CaloShowerROMaker] found energy (energy corr) (edep_npe) / nStepsMC / nPE "
                                   <<diagSum.totEdep<<" ("<<diagSum.totEdepCorr<<") ("<<diagSum.totEdepNPE<<") / "<<diagSum.totSteps
                                   <<" / "<<diagSum.totNPE<<std::endl;
  }


  //----------------------------------------------------------------------------------------------------------------------------------
  // Energy = (1+a*(Z/L-1/2))*energy where Z is position along the crystal, L the crystal length,
  // and a the non-uniformity level (response(z=1)-response(z=0))/response(z=1/2).
  float CaloShowerROMaker::LRUCorrection(float normalizedPosZ, float edepInit, float lru) const
  {
      const float edep = edepInit * (lru*(normalizedPosZ - 0.5f) + 1.0f);
      if (diagLevel_ > 2) std::cout<<"[CaloShowerROMaker::LRUCorrection] before / after LRU -> edep_corr = "
                                   << edepInit<<"  /  "<<edep<<"  at position Z="<<normalizedPosZ<<std::endl;
      return edep;
  }

  //-------------------------------------------------------------------------------------------------------------------------------------------------------
  void CaloShowerROMaker::dumpCaloShowerSim(const CaloShowerSimCollection& caloShowerSims) const
  {
       std::cout<<"[CaloShowerROMaker] Checking Sims"<<std::endl;
       float csmEtot(0);
       for (const auto& csm : caloShowerSims)
       {
           csmEtot += csm.energyDep();
           std::cout<<csm.crystalID()<<" "<<csm.sim()<<" "<<csm.time()<<" "<<csm.energyDep()<<" "<<csm.energyDepG4()<<std::endl;
           for (const auto& st : csm.caloShowerSteps()) std::cout<<"  "<<st<<std::endl;
       }
       std::cout<<"[CaloShowerROMaker] CSM Etot "<<csmEtot<<std::endl;
  }

}
DEFINE_ART_MODULE(mu2e::CaloShowerROMaker)
