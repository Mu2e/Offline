//
// Create a compressed representation of Calorimeter StepPointMCs.
//
// We recollect all StepPointMCs attached to a SimParticle ancestor and collapse them into
// CaloShowerStep objects. A SimParticle entering the calorimeter is an "ancestor" SimParticle;
// it generates a shower of SimParticles and StepPointMCs in the crystal, which is compressed.
// At the end all StepPointMCs can be dropped, as well as a large fraction of SimParticles, with
// negligible loss of information.
//
// Note: if a SimParticle enters the calorimeter and generates a secondary that hits another disk
// (e.g. a photon leaks from disk 1 into disk 2), the SimParticle hitting the second disk is itself
// considered an ancestor.
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "Offline/CalorimeterGeom/inc/Calorimeter.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/CaloMC/inc/ShowerStepUtil.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "cetlib_except/exception.h"

#include "CLHEP/Vector/ThreeVector.h"

#include <algorithm>
#include <iostream>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>


namespace {

   // Accumulates the steps of one ancestor SimParticle together with the set of
   // (secondary) SimParticles that were folded into it.
   class CaloCompressUtil
   {
       public:
           using SimPtr = art::Ptr<mu2e::SimParticle>;

           const std::vector<const mu2e::StepPointMC*>& steps() const {return steps_;}
           const std::set<SimPtr>&                      sims()  const {return sims_;}

           void fill(const mu2e::StepPointMC* step, const std::vector<SimPtr>& sims)
           {
               steps_.push_back(step);
               sims_.insert(sims.begin(), sims.end());
           }

        private:
           std::vector<const mu2e::StepPointMC*> steps_;
           std::set<SimPtr>                      sims_;
   };

   struct DiagSummary
   {
       float    totalEdep_{0.f};
       unsigned totalStep_{0}, totalSim_{0}, totalChk_{0}, nCompress_{0}, nCompressInfo_{0};
       void reset() {*this = DiagSummary{};}
   };

}


namespace mu2e {

  class CaloShowerStepMaker : public art::EDProducer
  {
      public:
         struct Config
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             fhicl::Sequence<std::string> caloStepPointCollection { Name("caloStepPointCollection"), Comment("Calo crystal stepPointMC collection name") };
             fhicl::Atom<unsigned>        numZSlices              { Name("numZSlices"),              Comment("Number of crystal longitudinal slices") };
             fhicl::Atom<float>           deltaTime               { Name("deltaTime"),               Comment("Max time difference to be inside a ShowerStep") };
             fhicl::Atom<bool>            compressData            { Name("compressData"),            Comment("Compress stepPointMC and SimParticles in crystal") };
             fhicl::Atom<double>          eDepThreshold           { Name("eDepThreshold"),           Comment("Threshold on energy deposited by SimParticle to keep it") };
             fhicl::Atom<int>             diagLevel               { Name("diagLevel"),               Comment("Debug"),0 };
         };

         explicit CaloShowerStepMaker(const art::EDProducer::Table<Config>& config);
         void produce(art::Event& e) override;


      private:
         using SimPtr      = art::Ptr<SimParticle>;
         using StepHandles = std::vector<art::ValidHandle<StepPointMCCollection>>;

         void makeCompressedHits      (const StepHandles&, CaloShowerStepCollection&, SimParticlePtrCollection&);
         void collectStepBySimAncestor(const Calorimeter&, const StepHandles&, std::map<SimPtr,CaloCompressUtil>&);
         void compressSteps           (const Calorimeter&, CaloShowerStepCollection&, int volId, const SimPtr&, std::vector<const StepPointMC*>& steps);
         void dumpAllInfo             (const StepHandles&, const Calorimeter&) const;

         std::vector<art::ProductToken<StepPointMCCollection>> stepTokens_;
         unsigned    numZSlices_;
         double      deltaTime_;
         bool        compressData_;
         double      eDepThreshold_;
         int         diagLevel_;
         double      zSliceSize_{0};
         DiagSummary diag_;
  };



  CaloShowerStepMaker::CaloShowerStepMaker(const art::EDProducer::Table<Config>& config) :
     art::EDProducer{config},
     numZSlices_   (config().numZSlices()),
     deltaTime_    (config().deltaTime()),
     compressData_ (config().compressData()),
     eDepThreshold_(config().eDepThreshold()),
     diagLevel_    (config().diagLevel())
  {
      if (numZSlices_ == 0)
         throw cet::exception("CONFIG") << "[CaloShowerStepMaker] numZSlices must be > 0\n";

      for (const auto& tag : config().caloStepPointCollection())
         stepTokens_.push_back(consumes<StepPointMCCollection>(art::InputTag(tag)));

      produces<CaloShowerStepCollection>();
      produces<SimParticlePtrCollection>();
  }


  //------------------------------------------------------------------------------------------------------------
  void CaloShowerStepMaker::produce(art::Event& event)
  {
      diag_.reset();
      if (diagLevel_ > 0) std::cout << "[CaloShowerStepMaker::produce] begin" << std::endl;

      auto caloShowerStepMCs = std::make_unique<CaloShowerStepCollection>();
      auto simsToKeep        = std::make_unique<SimParticlePtrCollection>();

      StepHandles crystalStepsHandles;
      crystalStepsHandles.reserve(stepTokens_.size());
      for (const auto& token : stepTokens_) crystalStepsHandles.push_back(event.getValidHandle(token));

      makeCompressedHits(crystalStepsHandles, *caloShowerStepMCs, *simsToKeep);

      event.put(std::move(caloShowerStepMCs));
      event.put(std::move(simsToKeep));

      if (diagLevel_ > 0) std::cout << "[CaloShowerStepMaker::produce] end" << std::endl;
  }


  //------------------------------------------------------------------------------------------------------------------
  void CaloShowerStepMaker::makeCompressedHits(const StepHandles& crystalStepsHandles,
                                               CaloShowerStepCollection& caloShowerStepMCs, SimParticlePtrCollection& simsToKeep)
  {
      const Calorimeter& cal = *(GeomHandle<Calorimeter>());
      zSliceSize_            = cal.G4Info().get<double>("crystalZLength")/float(numZSlices_) + 1e-5;

      // Collect the StepPointMCs produced by each ancestor SimParticle
      std::map<SimPtr,CaloCompressUtil> crystalAncestorsMap;
      collectStepBySimAncestor(cal, crystalStepsHandles, crystalAncestorsMap);

      if (diagLevel_ > 2) dumpAllInfo(crystalStepsHandles, cal);

      // Loop over ancestor SimParticles, check if they are compressible, and produce the corresponding CaloShowerStep
      std::set<SimPtr> simsToKeepUnique;
      for (const auto& [sim, info] : crystalAncestorsMap)
      {
          diag_.totalSim_ += info.sims().size();

          std::map<unsigned,std::vector<const StepPointMC*>> crystalMap;
          for (const StepPointMC* step : info.steps()) crystalMap[step->volumeId()].push_back(step);

          for (auto& [crid, steps] : crystalMap)
          {
              //Filter very small energy deposits at this stage
              double eDep(0);
              for (const auto* step : steps) eDep += step->totalEDep();
              if (eDep < eDepThreshold_) continue;

              if (compressData_)
              {
                  simsToKeepUnique.insert(sim);
                  compressSteps(cal, caloShowerStepMCs, crid, sim, steps);
              }
              else
              {
                  std::map<SimPtr, std::vector<const StepPointMC*>> newSimStepMap;
                  for (const StepPointMC* step : steps) newSimStepMap[step->simParticle()].push_back(step);
                  for (auto& [stepSim, stepVec] : newSimStepMap)
                  {
                      compressSteps(cal, caloShowerStepMCs, crid, stepSim, stepVec);
                      simsToKeepUnique.insert(stepSim);
                  }
              }
          }
          ++diag_.nCompressInfo_;
          if (compressData_) ++diag_.nCompress_;
      }

      simsToKeep.assign(simsToKeepUnique.begin(), simsToKeepUnique.end());

      // Final diag info
      if (diagLevel_ > 1)
      {
          std::cout<<"CaloShowerStepMaker summary"<<std::endl;

          std::set<int> volIds;
          for (const auto& css : caloShowerStepMCs) volIds.insert(css.volumeG4ID());

          for (int volId : volIds)
          {
             std::map<art::Ptr<SimParticle>, double> simMap;
             for (const auto& css : caloShowerStepMCs)
                if (css.volumeG4ID()==volId) simMap[css.simParticle()] += css.energyDepG4();

             for (const auto& [simPtr, energy] : simMap)
                std::cout<<"Vol id: "<<volId<<"  Sim id: "<<simPtr.id()<<"   energy="<<energy<<std::endl;
          }
      }

      if (diagLevel_ > 0)
        std::cout << "[CaloShowerStepMaker] compressed "<<diag_.nCompress_<<" / "<<diag_.nCompressInfo_<<" incoming SimParticles"<<std::endl
                  << "[CaloShowerStepMaker] keeping "<<simsToKeep.size()<<" SimParticles"<<std::endl
                  << "[CaloShowerStepMaker] Total sims init: "<<diag_.totalSim_<<std::endl
                  << "[CaloShowerStepMaker] Total caloShower steps: "<<caloShowerStepMCs.size()<<std::endl
                  << "[CaloShowerStepMaker] Total energy deposited / number of stepPointMC: "<<diag_.totalEdep_<<" / "<<diag_.totalStep_<<std::endl
                  << "[CaloShowerStepMaker] Total stepPointMCs seen: "<<diag_.totalChk_<<std::endl;
  }


  //------------------------------------------------------------------------------------------------------------------
  void CaloShowerStepMaker::collectStepBySimAncestor(const Calorimeter& cal,
                                                     const StepHandles& stepsHandles,
                                                     std::map<SimPtr,CaloCompressUtil>& ancestorsMap)
  {
     std::vector<SimPtr>               inspectedSims;
     std::unordered_map<SimPtr,SimPtr> simToAncestorMap;

     for (const auto& handle : stepsHandles)
     {
         const StepPointMCCollection& steps(*handle);

         for (const auto& step : steps)
         {
             SimPtr sim = step.simParticle();

             inspectedSims.clear();
             while (sim->hasParent())
             {
                 const auto alreadyInspected = simToAncestorMap.find(sim);
                 if (alreadyInspected != simToAncestorMap.end()) {sim = alreadyInspected->second; break;}
                 inspectedSims.push_back(sim);

                 if (!cal.isInsideAnyCrystal(sim->startPosition()))                    break;
                 if (!cal.isInsideSameDisk(sim->startPosition(),sim->endPosition()))   break;

                 sim = sim->parent();
             }

             for (const SimPtr& inspectedSim : inspectedSims) simToAncestorMap[inspectedSim] = sim;
             ancestorsMap[sim].fill(&step, inspectedSims);

             diag_.totalEdep_ += step.totalEDep();
         }
         diag_.totalStep_ += steps.size();
      }
  }


  //-------------------------------------------------------------------------------------------------------------------------------
  void CaloShowerStepMaker::compressSteps(const Calorimeter& cal, CaloShowerStepCollection& caloShowerStepMCs,
                                          int volId, const SimPtr& sim, std::vector<const StepPointMC*>& steps)
  {
     std::sort(steps.begin(), steps.end(), [](const StepPointMC* a, const StepPointMC* b){return a->time() < b->time();});

     ShowerStepUtil buffer(numZSlices_, ShowerStepUtil::weight_type::energy);

     for (const StepPointMC* step : steps)
     {
         const CLHEP::Hep3Vector pos = cal.mu2eToCrystal(volId, step->position());

         // clamp: a step at/beyond the back face would otherwise index past the last slice
         const unsigned idx = std::min<unsigned>(unsigned(std::max(1e-6, pos.z())/zSliceSize_), numZSlices_-1);

         if (buffer.entries(idx)>0 && (step->time()-buffer.t0(idx) > deltaTime_))
         {
             if (diagLevel_ > 2) {std::cout<<"[CaloShowerStepMaker::compressSteps] inserted  "; buffer.printBucket(idx);}
             diag_.totalChk_ += buffer.entries(idx);

             caloShowerStepMCs.emplace_back(volId, sim, buffer.entries(idx), buffer.time(idx), buffer.energyG4(idx),
                                            buffer.energyVis(idx), buffer.pIn(idx), buffer.pos(idx));
             buffer.reset(idx);
         }

         buffer.add(idx, step->totalEDep(), step->visibleEDep(), step->time(), step->momentum().mag(), pos);
     }

     //do not forget to flush the final buffer(s) :-)
     for (unsigned i=0; i<buffer.nBuckets(); ++i)
     {
         if (buffer.entries(i) == 0) continue;

         if (diagLevel_ > 2) {std::cout<<"[CaloShowerStepMaker::compressSteps] inserted "; buffer.printBucket(i);}
         diag_.totalChk_ += buffer.entries(i);

         caloShowerStepMCs.emplace_back(volId, sim, buffer.entries(i), buffer.time(i), buffer.energyG4(i),
                                        buffer.energyVis(i), buffer.pIn(i), buffer.pos(i));
     }
  }

  //-------------------------------------------------------------------------------------------------------------
  void CaloShowerStepMaker::dumpAllInfo(const StepHandles& stepsHandles, const Calorimeter& cal) const
  {
      std::cout<<"Dumping StepPointMCs  Mu2e / crystal / disk / diskFF frames"<<std::endl;
      for (const auto& handle : stepsHandles)
      {
          const StepPointMCCollection& steps(*handle);
          std::cout<<steps.size()<<std::endl;
          for (const auto& step : steps)
            std::cout<<step.volumeId()<<" "<<step.totalEDep()<<" "<<step.position()<<" "
                     <<cal.mu2eToCrystal(step.volumeId(),step.position())<<"   "
                     <<cal.mu2eToDisk(cal.crystal(step.volumeId()).diskID(),step.position())<<"   "
                     <<cal.mu2eToDiskFF(cal.crystal(step.volumeId()).diskID(),step.position())<<std::endl;
      }
  }

}

DEFINE_ART_MODULE(mu2e::CaloShowerStepMaker)
