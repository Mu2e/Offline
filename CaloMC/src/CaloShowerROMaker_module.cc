//
// Transform the energy deposited in the scintillator into photo-electrons (PE) seen by the photosensor. 
// Includes corrections from Birks law, longitudinal response uniformity and photo-statistcs fluctuations.
// The PE are generated individually and corrected for transit time.
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"

#include "CaloMC/inc/CaloPhotonPropagation.hh"
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/CaloShowerStep.hh"
#include "MCDataProducts/inc/CaloShowerRO.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "SeedService/inc/SeedService.hh"
#include "CLHEP/Random/RandPoissonQ.h"
#include "CLHEP/Random/RandFlat.h"

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <utility>

#include "TFile.h"
#include "TH2.h"

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
         timeCorr_ = std::min(timeCorr,timeCorr_);
      }

      std::vector<art::Ptr<mu2e::CaloShowerStep>> steps_;
      float edepCorr_,timeCorr_;
  };
  
  struct diagSummary
  {
      diagSummary() : totSteps(0),totNPE(0),totEdep(0.),totEdepCorr(0.),totEdepNPE(0.) {};
      int   totSteps,totNPE;
      float totEdep,totEdepCorr,totEdepNPE;
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
             using SPTO    = SimParticleTimeOffset::Config;
             fhicl::Sequence<art::InputTag>  caloShowerStepCollection { Name("caloShowerStepCollection"), Comment("Compressed shower inputs for calo crystals") };
             fhicl::Table<SPTO>              timeOffsets              { Name("TimeOffsets"),              Comment("Time maps to apply to sim particles before digitization.") };
             fhicl::Atom<float>              blindTime                { Name("blindTime"),                Comment("Minimum time of hit to be digitized") };
             fhicl::Atom<bool>               LRUCorrection            { Name("LRUCorrection"),            Comment("Include LRU corrections") };
             fhicl::Atom<bool>               BirksCorrection          { Name("BirksCorrection"),          Comment("Include Birks corrections") };
             fhicl::Atom<bool>               PEStatCorrection         { Name("PEStatCorrection"),         Comment("Include PE Poisson fluctuations") };
             fhicl::Atom<bool>               addTravelTime            { Name("addTravelTime"),            Comment("Include light propagation time") };
             fhicl::Atom<int>                diagLevel                { Name("diagLevel"),                Comment("Diag Level"),0 };
         };

         explicit CaloShowerROMaker(const art::EDProducer::Table<Config>& config) :
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
            photonProp_       (engine_)
         {
             // the following consumes statements are necessary because SimParticleTimeOffset::updateMap calls getValidHandle.
             for (auto const& tag : config().caloShowerStepCollection()) crystalShowerTokens_.push_back(consumes<CaloShowerStepCollection>(tag));
             for (auto const& tag : config().timeOffsets().inputs()) consumes<SimParticleTimeMap>(tag);
             produces<CaloShowerROCollection>();
             produces<CaloShowerSimCollection>();
         }

         void beginJob() override;
         void beginRun(art::Run& aRun) override;
         void produce(art::Event& e) override;


      private:
         using StepHandles = std::vector<art::ValidHandle<CaloShowerStepCollection>>;

         void  makeReadoutHits   (const StepHandles&, CaloShowerROCollection&, CaloShowerSimCollection&);
         float LRUCorrection     (int, float, float, const ConditionsHandle<CalorimeterCalibrations>&);
         float PECorrection      (int, float, float);
         void  dumpCaloShowerSim (const CaloShowerSimCollection& caloShowerSims);

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
         CaloPhotonPropagation   photonProp_;
         TH2F*                   hTime_;
         TH1F*                   hEtot_;
         TH1F*                   hECorrtot_;
         TH1F*                   hStot_;

  };


  //-----------------------------------------------
  void CaloShowerROMaker::beginJob()
  {      
      if (diagLevel_ > 1)
      {
          art::ServiceHandle<art::TFileService> tfs;
          hTime_      = tfs->make<TH2F>("hTIme",     "Photon prop time",       40,0,200,50,0,20);
          hEtot_      = tfs->make<TH1F>("hEtot",     "Total E dep",           150,     0,   150);
          hECorrtot_  = tfs->make<TH1F>("hECorrtot", "Total E dep corr",      150,     0,   150);
          hStot_      = tfs->make<TH1F>("hStot",     "Total Compress steps",  100,     0,   1000);
      }
  }
  
  //-----------------------------------------------
  void CaloShowerROMaker::beginRun(art::Run& aRun)
  {      
      photonProp_.buildTable();
  }


  //---------------------------------------------------------------
  void CaloShowerROMaker::produce(art::Event& event)
  {
      if (diagLevel_ > 0) std::cout << "[CaloShowerROMaker::produce] begin" << std::endl;

      //update condition cache
      ConditionsHandle<AcceleratorParams> accPar("ignored");
      mbtime_ = accPar->deBuncherPeriod;
      toff_.updateMap(event);

      // Containers to hold the output hits.
      auto CaloShowerROs  = std::make_unique<CaloShowerROCollection>();
      auto caloShowerSims = std::make_unique<CaloShowerSimCollection>();

      StepHandles newCrystalShowerTokens;
      std::transform(std::begin(crystalShowerTokens_), std::end(crystalShowerTokens_),  
                     back_inserter(newCrystalShowerTokens), 
                     [&event](const auto& token) {return event.getValidHandle(token);});
      
      makeReadoutHits(newCrystalShowerTokens, *CaloShowerROs, *caloShowerSims);

      // Add the output hit collection to the event
      event.put(std::move(CaloShowerROs));
      event.put(std::move(caloShowerSims));

      if (diagLevel_ > 0) std::cout << "[CaloShowerROMaker::produce] end" << std::endl;
  }



  //-----------------------------------------------------------------------------------------------------
  void CaloShowerROMaker::makeReadoutHits(const StepHandles& crystalShowerHandles,
                                          CaloShowerROCollection& CaloShowerROs, 
                                          CaloShowerSimCollection& caloShowerSims)
  {

      GlobalConstantsHandle<ParticleDataTable>  pdt;
      ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

      const Calorimeter& cal       = *(GeomHandle<Calorimeter>());
      const int   nROs             = cal.caloInfo().getInt("nSiPMPerCrystal");
      const float cryhalflength    = cal.caloInfo().getDouble("crystalZLength")/2.0;

      std::map<int,std::vector<StepEntry>> simEntriesMap;

      //-----------------------------------------------------------------------
      //store corrected energy deposits for each redouts
      
      diagSummary diagSum;

      for (const auto& showerHandle: crystalShowerHandles)
      {
          const CaloShowerStepCollection& caloShowerSteps(*showerHandle);
          for (auto istep = caloShowerSteps.begin(); istep !=caloShowerSteps.end(); ++istep)
          {
              const CaloShowerStep& step = *istep;

              // time folding and filtering, see docdb-3425 for a stunning explanation
              double hitTimeUnfolded = toff_.totalTimeOffset(istep->simParticle())+istep->time();
              double hitTime         = fmod(hitTimeUnfolded,mbtime_);

              if (hitTime < blindTime_ || hitTime > mbtime_ ) continue;

              size_t idx = std::distance(caloShowerSteps.begin(), istep);
              art::Ptr<CaloShowerStep> stepPtr = art::Ptr<CaloShowerStep>(showerHandle,idx);            
              
              int   crystalID  = step.volumeG4ID();
              int   SiPMIDBase = cal.caloIDMapper().SiPMIDFromCrystalID(crystalID);              
              float posZ       = step.position().z();
              
              float edep_corr(step.energyDepG4());
              if (BirksCorrection_) edep_corr = step.energyDepBirks();
              if (LRUCorrection_)   edep_corr = LRUCorrection(crystalID, posZ/cryhalflength, edep_corr, calorimeterCalibrations);

              // Generate individual PEs and their arrival times
              for (int i=0; i<nROs; ++i)
              {
                  int SiPMID = SiPMIDBase + i;
                  float peMeV = calorimeterCalibrations->peMeV(SiPMID);
                  int NPE     = randPoisson_.fire(edep_corr*peMeV);                  
                  if (NPE==0) continue;
                  
                  std::vector<float> PETime(NPE,hitTime);                  
                  if (addTravelTime_) 
                  {
                      for (auto& time : PETime) time += photonProp_.propTimeSimu(2.0*cryhalflength-posZ);
                  }    
                  CaloShowerROs.push_back(CaloShowerRO(SiPMID,stepPtr,PETime));                  
                                    
                  if (diagLevel_ > 2) std::cout<<"[CaloShowerROMaker::generatePE] SiPMID:"<<SiPMID<<"  energy / NPE = "<<edep_corr<<"  /  "<<NPE<<std::endl;
                  if (diagLevel_ > 2) {std::cout<<"Time hit "<<std::endl; for (auto time : PETime) std::cout<<time<<" "; std::cout<<std::endl;}
                  if (diagLevel_ > 1) for (const auto& time : PETime) hTime_->Fill(2.0*cryhalflength-posZ,time-hitTime);
                                    
                  diagSum.totNPE     += NPE;
                  diagSum.totEdepNPE += double(NPE)/peMeV/2.0; //average between the two RO
              }
              diagSum.totEdep     += step.energyDepG4();
              diagSum.totEdepCorr += edep_corr;
              diagSum.totSteps    += step.nCompress();

              //Produce an MC object that include the step and additional information for each original step
              simEntriesMap[crystalID].push_back(StepEntry(stepPtr,edep_corr,hitTime));
          }
          
          auto sortFunctor = [](const auto& a, const auto& b){return a.SiPMID() < b.SiPMID();};
          std::sort(CaloShowerROs.begin(),CaloShowerROs.end(),sortFunctor); 
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
                 summaryMap.insert(std::make_pair(sim,SimParticleSummary(newStep.step_,newStep.edepCorr_,newStep.timeCorr_)));
              else
                 mfind->second.add(newStep.step_,newStep.edepCorr_,newStep.timeCorr_);
          }

          // create the CaloShowerSim (MC truth) objects for a given crystalID
          for (auto& kvsumm : summaryMap) caloShowerSims.push_back(CaloShowerSim(kvsumm.second.steps_, kvsumm.second.edepCorr_, kvsumm.second.timeCorr_));
      }
       


      //--------------------------------------------------
      // Diag
      if (diagLevel_ > 2) dumpCaloShowerSim(caloShowerSims);

      if (diagLevel_ > 1)
      {
         hEtot_->Fill(diagSum.totEdep);
         hECorrtot_->Fill(diagSum.totEdepCorr);
         hStot_->Fill(diagSum.totSteps);
         
         std::set<int> crIds;
         for (const auto& css : caloShowerSims) crIds.insert(css.crystalID());

         for (auto crId : crIds)
         {
            std::map<const art::Ptr<SimParticle>, double> simMap;
            for (const auto& css : caloShowerSims) if (css.crystalID()==crId) simMap[css.sim()] += css.energyDep();
            for (auto& kv : simMap) std::cout<<"CrId: "<<crId<<"  Sim id: "<<kv.first.id()<<"   energy="<<kv.second<<std::endl;   
         }
      }

      if (diagLevel_ > 0) std::cout<<"[CaloShowerROMaker::produce] found energy (energy corr) (edep_npe) / nStepsMC / nPE "
                                   <<diagSum.totEdep<<" ("<<diagSum.totEdepCorr<<") ("<<diagSum.totEdepNPE<<") / "<<diagSum.totSteps
                                   <<" / "<<diagSum.totNPE<<std::endl;

  }


  //----------------------------------------------------------------------------------------------------------------------------------
  // apply a correction of type Energy = ((1-s)*Z/HL+s)*energy where Z position along the crystal, HL is the crystal half-length
  // and s is the intercept at Z=0 (i.e. non-uniformity factor, e.g. 5% -> s = 1.05)
  float CaloShowerROMaker::LRUCorrection(int crystalID, float normalizedPosZ, float edepInit, 
                                         const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
      float alpha  = calorimeterCalibrations->LRUpar0(crystalID);
      float factor = (1.0-alpha)*normalizedPosZ +alpha;
      float edep   = edepInit*factor;

      if (diagLevel_ > 2) std::cout<<"[CaloShowerROMaker::LRUCorrection] before / after LRU -> edep_corr = "
                                   << edepInit<<"  /  "<<edep<<"  at position Z="<<normalizedPosZ<<std::endl;
      return edep;
  }

  //-------------------------------------------------------------------------------------------------------------------------------------------------------
  void CaloShowerROMaker::dumpCaloShowerSim(const CaloShowerSimCollection& caloShowerSims)
  {
       std::cout<<"[CaloShowerROMaker] Checking Sims"<<std::endl;
       float csmEtot(0);
       for (auto& csm :  caloShowerSims)
       {
           csmEtot += csm.energyDep();
           std::cout<<csm.crystalID()<<" "<<csm.sim()<<" "<<csm.time()<<" "<<csm.energyDep()<<" "<<csm.energyDepG4()<<std::endl;
           for (auto& st : csm.caloShowerSteps()) std::cout<<"  "<<st<<std::endl;
       }
       std::cout<<"[CaloShowerROMaker] CSM Etot "<<csmEtot<<std::endl;
  }

}

using mu2e::CaloShowerROMaker;
DEFINE_ART_MODULE(CaloShowerROMaker);
