//
// Transform the energy deposited in the scintillator into photo-electrons (PE) seen by the photosensor. 
// Includes corrections from Birks law, longitudinal response uniformity and photo-statistcs fluctuations.
// The PE are generated individually and corrected for transit time.
//
// Original author Bertrand Echenard
//
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepROCollection.hh"
#include "MCDataProducts/inc/CaloShowerSimCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "SeedService/inc/SeedService.hh"

#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandExponential.h"

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <utility>


namespace {

  struct StepEntry
  {
      StepEntry(const art::Ptr<mu2e::CaloShowerStep>& step, float edepCorr, float time) :
        step_(step), edepCorr_(edepCorr), time_(time) {}

      art::Ptr<mu2e::CaloShowerStep> step_;
      float edepCorr_;
      float time_;
  };

  struct SimParticleSummary
  {
      SimParticleSummary(const art::Ptr<mu2e::CaloShowerStep>& step, float edepCorr, float time) :
        steps_(), edepCorr_(edepCorr), time_(time) {steps_.push_back(step);}

      void add(const art::Ptr<mu2e::CaloShowerStep>& step, float edepCorr, float time)
      {
         steps_.push_back(step);
         edepCorr_ += edepCorr;
         time_ = std::min(time_,time);
      }

      std::vector<art::Ptr<mu2e::CaloShowerStep>> steps_;
      float edepCorr_;
      float time_;
  };
}




namespace mu2e {

  class CaloShowerStepROFromShowerStep : public art::EDProducer 
  {
      public:
         struct Config 
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             using SPTO    = SimParticleTimeOffset::Config;

             fhicl::Sequence<art::InputTag>  caloShowerStepCollection { Name("caloShowerStepCollection"), Comment("Compressed shower inputs for calo crystals") };
             fhicl::Table<SPTO>              timeOffsets              { Name("TimeOffsets"),             Comment("Time maps to apply to sim particles before digitization.") };
             fhicl::Atom<float>              blindTime                { Name("blindTime"),               Comment("Minimum time of hit to be digitized") };
             fhicl::Atom<bool>               LRUCorrection            { Name("LRUCorrection"),           Comment("Include LRU corrections") };
             fhicl::Atom<bool>               BirksCorrection          { Name("BirksCorrection"),         Comment("Include Birks corrections") };
             fhicl::Atom<bool>               PEStatCorrection         { Name("PEStatCorrection"),        Comment("Include PE Poisson fluctuations") };
             fhicl::Atom<bool>               addTravelTime            { Name("addTravelTime"),           Comment("Include light propagation time") };
             fhicl::Atom<int>                diagLevel                { Name("diagLevel"),               Comment("Diag Level"),0 };
         };

         explicit CaloShowerStepROFromShowerStep(const art::EDProducer::Table<Config>& config) :
            EDProducer{config},
            toff_             (config().timeOffsets()),
            blindTime_        (config().blindTime()),
            LRUCorrection_    (config().LRUCorrection()),
            BirksCorrection_  (config().BirksCorrection()),
            PEStatCorrection_ (config().PEStatCorrection()),
            addTravelTime_    (config().addTravelTime()),
            diagLevel_        (config().diagLevel()),
            engine_           (createEngine(art::ServiceHandle<SeedService>()->getSeed())),
            randPoisson_      (engine_),
            randGauss_        (engine_),
            randExpo_         (engine_)
         {
             // the following consumes statements are necessary because SimParticleTimeOffset::updateMap calls getValidHandle.
             for (auto const& tag : config().caloShowerStepCollection()) crystalShowerTokens_.push_back(consumes<CaloShowerStepCollection>(tag));
             for (auto const& tag : config().timeOffsets().inputs()) consumes<SimParticleTimeMap>(tag);
             produces<CaloShowerStepROCollection>();
             produces<CaloShowerSimCollection>();
         }

         void beginJob() override;
         void produce(art::Event& e) override;


      private:
         using StepHandles = std::vector<art::ValidHandle<CaloShowerStepCollection>>;

         void  makeReadoutHits  (const StepHandles&, CaloShowerStepROCollection&, CaloShowerSimCollection&);
         float LRUCorrection    (int crystalID, float normalizedPosZ, float edepInit, const ConditionsHandle<CalorimeterCalibrations>&);
         float PECorrection     (int crystalID, float edepInit, float NpePerMeV);
         void  dumpCaloShowerSim(const CaloShowerSimCollection& caloShowerSims);

         std::vector<art::ProductToken<CaloShowerStepCollection>> crystalShowerTokens_;
         SimParticleTimeOffset   toff_;
         float                   blindTime_;
         float                   mbtime_;
         bool                    LRUCorrection_;
         bool                    BirksCorrection_;
         bool                    PEStatCorrection_;
         bool                    addTravelTime_;
         int                     diagLevel_;
         CLHEP::HepRandomEngine& engine_;
         CLHEP::RandPoissonQ     randPoisson_;
         CLHEP::RandGaussQ       randGauss_;
         CLHEP::RandExponential  randExpo_;

  };


  //-----------------------------------------------
  void CaloShowerStepROFromShowerStep::beginJob()
  {}


  //---------------------------------------------------------------
  void CaloShowerStepROFromShowerStep::produce(art::Event& event)
  {
      if (diagLevel_ > 0) std::cout << "[CaloShowerStepROFromShowerStep::produce] begin" << std::endl;

      //update condition cache
      ConditionsHandle<AcceleratorParams> accPar("ignored");
      mbtime_ = accPar->deBuncherPeriod;
      toff_.updateMap(event);

      // Containers to hold the output hits.
      auto caloShowerStepROs = std::make_unique<CaloShowerStepROCollection>();
      auto caloShowerSims    = std::make_unique<CaloShowerSimCollection>();

      StepHandles newCrystalShowerTokens;
      std::transform(std::begin(crystalShowerTokens_), std::end(crystalShowerTokens_),  back_inserter(newCrystalShowerTokens), 
                     [&event](const auto& token) {return event.getValidHandle(token);});
      
      makeReadoutHits(newCrystalShowerTokens, *caloShowerStepROs, *caloShowerSims);

      // Add the output hit collection to the event
      event.put(std::move(caloShowerStepROs));
      event.put(std::move(caloShowerSims));

      if (diagLevel_ > 0) std::cout << "[CaloShowerStepROFromShowerStep::produce] end" << std::endl;
  }



  //-----------------------------------------------------------------------------------------------------
  void CaloShowerStepROFromShowerStep::makeReadoutHits(const StepHandles& crystalShowerHandles,
                                                       CaloShowerStepROCollection& caloShowerStepROs, 
                                                       CaloShowerSimCollection& caloShowerSims)
  {

      GlobalConstantsHandle<ParticleDataTable>  pdt;
      ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

      const Calorimeter& cal       = *(GeomHandle<Calorimeter>());
      const int   nROs             = cal.caloInfo().nROPerCrystal();
      const float cryhalflength    = cal.caloInfo().getDouble("crystalZLength")/2.0;
      const float refractiveIndex  = cal.caloInfo().getDouble("refractiveIndex");
      const float lightSpeed       = 300; // mm/ns
      //const float crystalDecayTime = cal.caloInfo().crystalDecayTime();

      std::map<int,std::vector<StepEntry>> simEntriesMap;
      std::map<int,float> totCrystalEnerMap;


      //-----------------------------------------------------------------------
      //store corrected energy deposits for each redouts

      int   totalSteps(0),totalNPE(0);
      float totalEdep(0.0),totalEdepCorr(0.0), totalEdepNPE(0.0);

      for (const auto& showerHandle: crystalShowerHandles)
      {
          const CaloShowerStepCollection& caloShowerSteps(*showerHandle);
          for (auto istep = caloShowerSteps.begin(); istep !=caloShowerSteps.end(); ++istep)
          {
              const CaloShowerStep& step = *istep;

              // time folding and filtering, see docdb-3425 for a stunning explanation
              float hitTimeUnfolded = toff_.totalTimeOffset(istep->simParticle())+istep->time();
              float hitTime         = fmod(hitTimeUnfolded,mbtime_);

              if (hitTime < blindTime_ || hitTime > mbtime_ ) continue;

              size_t idx = std::distance(caloShowerSteps.begin(), istep);
              art::Ptr<CaloShowerStep> stepPtr = art::Ptr<CaloShowerStep>(showerHandle,idx);            
              
              int   crystalID = step.volumeId();
              int   ROIDBase  = cal.caloInfo().ROBaseByCrystal(crystalID);              
              float posZ      = step.position().z();
              
              float edep_corr(step.energyDepG4());
              if (BirksCorrection_) edep_corr = step.energyDepBirks();
              if (LRUCorrection_)   edep_corr = LRUCorrection(crystalID, posZ/cryhalflength, edep_corr, calorimeterCalibrations);

              // Generate individual PEs and their arrival times
              for (int i=0; i<nROs; ++i)
              {
                  int ROID = ROIDBase + i;

                  int NPE = randPoisson_.fire(edep_corr*calorimeterCalibrations->peMeV(ROID));
                  if (diagLevel_ > 3) std::cout<<"[CaloShowerStepROFromShowerStep::generatePE] energy / NPE = "<<edep_corr<<"  /  "<<NPE<<"  "<<edep_corr*calorimeterCalibrations->peMeV(ROID)<<std::endl;
                  if (NPE==0) continue;
                  
                  std::vector<float> PETime(NPE,hitTime);
                  for (auto& time : PETime) time += (2.0*cryhalflength-posZ)*refractiveIndex/lightSpeed;
                  
                  caloShowerStepROs.push_back(CaloShowerStepRO(ROID,stepPtr,PETime));                  
                  
                  if (diagLevel_ > 3) {std::cout<<"Time hit "<<std::endl; for (auto time : PETime) std::cout<<time<<" "; std::cout<<std::endl;}
                  totalNPE += NPE;
                  totalEdepNPE +=NPE/calorimeterCalibrations->peMeV(ROID);
             }


              //Produce an object that include the step and additional information for each original Step
              simEntriesMap[crystalID].push_back(StepEntry(stepPtr,edep_corr,hitTime));

              if (diagLevel_ > 0)
              {
                  totCrystalEnerMap[crystalID] += edep_corr;
                  totalEdep     += step.energyDepG4();
                  totalEdepCorr += edep_corr;
                  totalSteps    += step.nCompress();
              }
          } 
      } 


      if (diagLevel_ > 1)
      {
          std::cout<<"[CaloShowerStepROFromShowerStep::produce] Crystal summary "<<std::endl;
          for (auto& kv : totCrystalEnerMap) std::cout<<"roId ="<< kv.first*2<<" "<<kv.second<<std::endl;
          if (diagLevel_ > 2) std::cout<<"[CaloShowerStepROFromShowerStep::produce] caloShowerStepROs summary"<<std::endl;
          if (diagLevel_ > 2) for (const auto& cst: caloShowerStepROs) std::cout<<cst.ROID()<<" "<<cst.NPE()<<std::endl;
      }



      //--------------------------------------------------
      // Produce the final MC truth info collecting energy deposits for each SimParticle in each crystal
      for (auto& kv : simEntriesMap)
      {
          std::vector<StepEntry> newSteps = kv.second;

          // fill the summary map for each simPtr for a given crystalID
          std::map<art::Ptr<SimParticle>,SimParticleSummary> summaryMap;
          for (auto& newStep : newSteps)
          {
              const art::Ptr<SimParticle>& sim = newStep.step_->simParticle();
              auto mfind = summaryMap.find(sim);
              if (mfind==summaryMap.end())
                 summaryMap.insert(std::make_pair(sim,SimParticleSummary(newStep.step_,newStep.edepCorr_,newStep.time_)));
              else
                 mfind->second.add(newStep.step_,newStep.edepCorr_,newStep.time_);
          }

          // create the CaloShowerSim (MC truth) objects for a given crystalID
          for (auto& kvsumm : summaryMap) caloShowerSims.push_back(CaloShowerSim(kvsumm.second.steps_, kvsumm.second.time_,kvsumm.second.edepCorr_));
      }


      if (diagLevel_ > 3) dumpCaloShowerSim(caloShowerSims);
      if (diagLevel_ > 0) std::cout<<"[CaloShowerStepROFromShowerStep::produce] found energy (energy corr) (edep_npe) / nStepsMC / nPE "
                                   <<totalEdep<<" ("<<totalEdepCorr<<") / ("<<totalEdepNPE<<") / "<<totalSteps<<" / "<<totalNPE<<std::endl;

  }


  //----------------------------------------------------------------------------------------------------------------------------------
  // apply a correction of type Energy = ((1-s)*Z/HL+s)*energy where Z position along the crystal, HL is the crystal half-length
  // and s is the intercept at Z=0 (i.e. non-uniformity factor, e.g. 5% -> s = 1.05)
  float CaloShowerStepROFromShowerStep::LRUCorrection(int crystalID, float normalizedPosZ, float edepInit, 
                                                      const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
      float alpha  = calorimeterCalibrations->LRUpar0(crystalID);
      float factor = (1.0-alpha)*normalizedPosZ +alpha;
      float edep   = edepInit*factor;

      if (diagLevel_ > 3) std::cout<<"[CaloShowerStepROFromShowerStep::LRUCorrection] before / after LRU -> edep_corr = "<< edepInit<<"  /  "<<edep<<"  at position Z="<<normalizedPosZ<<std::endl;
      return edep;
  }

  //-------------------------------------------------------------------------------------------------------------------------------------------------------
  void CaloShowerStepROFromShowerStep::dumpCaloShowerSim(const CaloShowerSimCollection& caloShowerSims)
  {
       std::cout<<"Checking Sims"<<std::endl;
       float csmEtot(0);
       for (auto& csm :  caloShowerSims)
       {
           csmEtot += csm.energyDep();
           std::cout<<csm.crystalId()<<" "<<csm.sim()<<" "<<csm.time()<<" "<<csm.energyDep()<<" "<<csm.energyDepG4()<<std::endl;
           for (auto& st : csm.caloShowerSteps()) std::cout<<"  "<<st<<std::endl;
       }
       std::cout<<"CSM Etot "<<csmEtot<<std::endl;
  }

}

using mu2e::CaloShowerStepROFromShowerStep;
DEFINE_ART_MODULE(CaloShowerStepROFromShowerStep);




/*
//little snipper I kept in case it will be needed...

float crystalDecayTime = cal.caloInfo().crystalDecayTime();

//std::vector<float> photonTime(nPhot);
//std::generate(photonTime.begin(),photonTime.end(),[&](){return hitTime+randExpo_.fire(crystalDecayTime);});
*/
