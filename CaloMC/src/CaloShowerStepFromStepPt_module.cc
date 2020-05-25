//
// An EDProducer Module that creates a compressed representation of Calorimeter StepPointMCs
//
// Original author Bertrand Echenard
//
// Basic idea: We recollect all StepPointMCS attached to a SimParticle ancestor. If the SimParticle is compressible, all StepPointsMC 
// are collapsed into CaloShowerStep objects. If not, we create a CaloShowerStep object for each SimParticle created by the "ancestor"
// SimParticle (basically, compress the StepPointMC for each SimParticle). We call a SimParticle entering the calorimeter an 
// Ancestore SimParticle. This ancestor will generate a shower of SimParticles and StepPointMcs in the crystal, which will be compressed
// At the end of the modules, all StepPointMCs can be dropped, as well as a large fraction of SimParticles 
// with almost no loss of information.
//
// Note: if a SimParticle enters the calorimeter, generates a secondary SimParticle that hit another section of the calorimeter 
// (e.g. a photon leaks from the first disk and hits the second disk), then the SimParticle hitting the second section is considered 
// to be an ancestor SimParticle
//
// The compressibility is determined by looking at the interaction codes of the StepPointMCs. These are currently hardcoded.
// Particles are compressed in small intervals of time and crystal longitudinal slices. There is an option to compress all particles.
//
//
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art_root_io/TFileService.h"
#include "art_root_io/TFileDirectory.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "CaloMC/inc/ShowerStepUtil.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Mu2eUtilities/inc/PhysicalVolumeMultiHelper.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "TH2F.h"

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <set>
#include <unordered_map>
#include <utility>



namespace {

   class CaloCompressUtil 
   {
       public:
           CaloCompressUtil() : steps_(), sims_(), procs_() {}

           const std::vector<const mu2e::StepPointMC*>& steps()        const {return steps_;}
           const std::set<art::Ptr<mu2e::SimParticle>>& sims()         const {return sims_;}
           const std::set<int>&                         processCodes() const {return procs_;}

           void fill(const mu2e::StepPointMC* step, std::vector<art::Ptr<mu2e::SimParticle>> sims)
           {
               steps_.push_back(step);
               procs_.insert(step->endProcessCode());               
               for (const auto& sim: sims) sims_.insert(sim); 
           }

        private:
           std::vector<const mu2e::StepPointMC*> steps_;
           std::set<art::Ptr<mu2e::SimParticle>> sims_;
           std::set<int>                         procs_;                              
   };


   struct diagSummary
   {
       diagSummary() : totalEdep_(0.0),totalStep_(0),totalSim_(0),totalChk_(0),nCompress_(0),nCompressAll_(0) {};
       void reset() {totalEdep_=0.0;totalStep_=totalSim_=totalChk_=nCompress_=nCompressAll_=0;}

       float     totalEdep_;
       unsigned  totalStep_,totalSim_,totalChk_,nCompress_,nCompressAll_;
   };

} 



namespace mu2e {

  class CaloShowerStepFromStepPt : public art::EDProducer 
  {
      public:      
         struct Config 
         {
             using Name    = fhicl::Name;
             using Comment = fhicl::Comment;
             fhicl::Atom<std::string>      caloStepPointCollection { Name("caloStepPointCollection"), Comment("Calo crystal stepPointMC collection name") };
             fhicl::Atom<art::InputTag>    physVolInfoInput        { Name("physVolInfoInput"),        Comment("Physics volume token names") };
             fhicl::Atom<unsigned>         numZSlices              { Name("numZSlices"),              Comment("Number of crystal longitudinal slices ") };
             fhicl::Atom<float>            deltaTime               { Name("deltaTime"),               Comment("Max time difference to be inside a ShowerStep") };
             fhicl::Atom<bool>             usePhysVolInfo          { Name("usePhysVolInfo"),          Comment("Use Physical Info volume names") };
             fhicl::Sequence<std::string>  caloMaterial            { Name("caloMaterial"),            Comment("List of calo material names") };
             fhicl::Atom<bool>             compressAll             { Name("compressAll"),             Comment("Compress all tracks") };
             fhicl::Atom<int>              diagLevel               { Name("diagLevel"),               Comment("Debug"),0 };
         };

         explicit CaloShowerStepFromStepPt(const art::EDProducer::Table<Config>& config);

         void beginJob() override;
         void beginSubRun(art::SubRun& sr) override;
         void produce( art::Event& e) override;


      private:
         using HandleVector = std::vector<art::Handle<StepPointMCCollection>>;
         using SimPtr       = art::Ptr<SimParticle>;
         
         void makeCompressedHits      (const HandleVector&, CaloShowerStepCollection&, SimParticlePtrCollection&);
         void collectStepBySimAncestor(const Calorimeter&, const PhysicalVolumeMultiHelper&, 
                                       const HandleVector&,std::map<SimPtr,CaloCompressUtil>&);
         void collectStepBySim        (const HandleVector&, std::map<SimPtr,std::vector<const StepPointMC*> >&);
         bool isInsideCalorimeter     (const Calorimeter& cal, const PhysicalVolumeMultiHelper&, const art::Ptr<SimParticle>&);
         bool isCompressible          (int, const std::set<int>&);
         void compressSteps           (const Calorimeter&, CaloShowerStepCollection&, int, const SimPtr&, std::vector<const StepPointMC*>&);
         void fillHisto1              (const Calorimeter&, const art::Ptr<SimParticle>&, const std::set<art::Ptr<SimParticle>>&);
         void fillHisto2              (int, float, const SimPtr&);
         void dumpAllInfo             (const HandleVector&, const Calorimeter&);


         
         std::string                                    calorimeterStepPoints_;
         art::InputTag                                  physVolInfoInput_;
         std::set<const PhysicalVolumeInfo*>  mapPhysVol_;
         bool                                           usePhysVol_;
         std::vector<std::string>                       caloMaterial_;
         int                                            numZSlices_;
         double                                         deltaTime_;
         bool                                           compressAll_;
         int                                            diagLevel_;
         std::map<int,std::set<int>>                    procCodes_;
         const PhysicalVolumeInfoMultiCollection*       vols_;
         double                                         zSliceSize_;

         diagSummary                                    diagSummary_;
         TH2F*                                          hStartPos_;
         TH2F*                                          hStopPos_;
         TH1F*                                          hStopPos2_;
         TH1F*                                          hStartPos2_;
         TH1F*                                          hZpos_;
         TH1F*                                          hEtot_;
         TH2F*                                          hZpos2_;
         TH1F*                                          hGenId_;
  };

  
  
  CaloShowerStepFromStepPt::CaloShowerStepFromStepPt(const art::EDProducer::Table<Config>& config) :
     art::EDProducer{config},     
     calorimeterStepPoints_(config().caloStepPointCollection()),
     physVolInfoInput_     (config().physVolInfoInput()),
     usePhysVol_           (config().usePhysVolInfo()),
     caloMaterial_         (config().caloMaterial()),
     numZSlices_           (config().numZSlices()),
     deltaTime_            (config().deltaTime()),
     compressAll_          (config().compressAll()), 
     diagLevel_            (config().diagLevel()),
     procCodes_(),
     vols_(),
     zSliceSize_(0),
     diagSummary_()
     {           
         consumesMany<StepPointMCCollection>();
         produces<CaloShowerStepCollection>();
         produces<SimParticlePtrCollection>();
     }


  //--------------------------------------------------------------------
  void CaloShowerStepFromStepPt::beginJob()
  {
      // procCodes are the process codes for the StepPointMC.
      // see MCDataProducts/inc/ProcessCode.hh for code numbering scheme
      //
      // --- These are hardcoded to make sure changes are intended and carefully considered ---
      //
      procCodes_[11].insert(   {2,12,13,16,17,21,23,29,40,49,58} );    // electron
      procCodes_[-11].insert(  {2,12,13,16,17,21,23,29,40,49,58} );    // positron
      procCodes_[22].insert(   {2,12,13,16,17,21,23,29,40,49,58} ); // photon
      procCodes_[2112].insert( {2,12,13,16,17,21,23,29,40,49,58,74} ); // neutron
      procCodes_[2212].insert( {16,17,21,23,29,40,45,49,58} );   // proton
      procCodes_[13].insert(   {2,12,13,16,17,21,23,29,30,31,34,40,49,58,59} ); // mu-
      procCodes_[-13].insert(  {2,12,13,16,17,21,23,29,30,31,34,40,49,58,59} ); // mu-

      if (diagLevel_ > 1)
      {
          art::ServiceHandle<art::TFileService> tfs;
          hStartPos_  = tfs->make<TH2F>("hStartPos", "Sim start position",  1000,  5000, 15000, 200, 0, 1000);
          hStopPos_   = tfs->make<TH2F>("hStopPos",  "Sim stop position",   1000,  5000, 15000, 200, 0, 1000);
          hStartPos2_ = tfs->make<TH1F>("hStartPos2","Sim stop position",   1000, 10000, 13000);
          hStopPos2_  = tfs->make<TH1F>("hStopPos2", "Sim stop position",   1000, 10000, 13000);
          hZpos_      = tfs->make<TH1F>("hZpos",     "Step z pos",            20,     0,    20);
          hZpos2_     = tfs->make<TH2F>("hZpos2",    "Step z pos",            20,     0,    20, 100, 0, 5);
          hEtot_      = tfs->make<TH1F>("hEtot",     "Sim stop position",    150,     0,   150);
          hGenId_     = tfs->make<TH1F>("hSimId",    "Gen Id",               150,    -10,  140);
      }
  }

  void CaloShowerStepFromStepPt::beginSubRun(art::SubRun& sr)
  {
      mapPhysVol_.clear();

      art::Handle<PhysicalVolumeInfoMultiCollection> volh;
      sr.getByLabel(physVolInfoInput_, volh);        
      if (!volh.isValid()) return;

      vols_ = volh.product();
      for (const auto& vol : *volh)
      {
          for (const auto& mv : vol.second) 
          {  
              if (std::find(caloMaterial_.begin(),caloMaterial_.end(), mv.second.materialName()) != caloMaterial_.end()) 
                 mapPhysVol_.insert(&mv.second);
          }
      }
  }


  //------------------------------------------------------------------------------------------------------------
  //
  void CaloShowerStepFromStepPt::produce(art::Event& event)
  {
      diagSummary_.reset();
      if (diagLevel_ > 0) std::cout << "[CaloShowerStepFromStepPt::produce] begin" << std::endl;

      //Create output collections
      auto caloShowerStepMCs = std::make_unique<CaloShowerStepCollection>();
      auto simsToKeep        = std::make_unique<SimParticlePtrCollection>();

      // Get the StepPointMCs from the event.
      art::ProductInstanceNameSelector getCrystalSteps(calorimeterStepPoints_);
      HandleVector crystalStepsHandles, readoutStepsHandles;
      event.getMany(getCrystalSteps, crystalStepsHandles);

      makeCompressedHits(crystalStepsHandles,*caloShowerStepMCs,*simsToKeep);

      // Add the output hit collection to the event
      event.put(std::move(caloShowerStepMCs));
      event.put(std::move(simsToKeep));

      if (diagLevel_ > 0) std::cout << "[CaloShowerStepFromStepPt::produce] end" << std::endl;
  }


  //------------------------------------------------------------------------------------------------------------------
  void CaloShowerStepFromStepPt::makeCompressedHits(const HandleVector& crystalStepsHandle, 
                                                    CaloShowerStepCollection& caloShowerStepMCs,SimParticlePtrCollection& simsToKeep)
  {
      PhysicalVolumeMultiHelper vi(*vols_);

      const Calorimeter& cal = *(GeomHandle<Calorimeter>());
      zSliceSize_            = cal.caloInfo().getDouble("crystalZLength")/float(numZSlices_)+1e-5;


      //-----------------------------------------------------------------
      // Collect the StepPointMC's produced by each SimParticle Ancestor
      std::map<SimPtr,CaloCompressUtil> crystalAncestorsMap;
      collectStepBySimAncestor(cal,vi,crystalStepsHandle,crystalAncestorsMap);
      
      if (diagLevel_ > 2) dumpAllInfo(crystalStepsHandle,cal);



      //---------------------------------------------------------------------------------------------------------------
      //Loop over ancestor simParticles, check if they are compressible, and produce the corresponding caloShowerStepMC

      std::set<SimPtr> SimsToKeepUnique;
      for (const auto& iter : crystalAncestorsMap )
      {
          const SimPtr&           sim  = iter.first;
          const CaloCompressUtil& info = iter.second;

          bool doCompress = isCompressible(sim->pdgId(),info.processCodes());
          diagSummary_.totalSim_ += info.sims().size();

          std::map<int,std::vector<const StepPointMC*>> crystalMap;
          for (const StepPointMC* step : info.steps()) crystalMap[step->volumeId()].push_back(step);

          for (const auto& iterCrystal : crystalMap)
          {
              int crid = iterCrystal.first;
              std::vector<const StepPointMC*> steps = iterCrystal.second;

              if (doCompress)
              {
                  SimsToKeepUnique.insert(sim);
                  compressSteps(cal, caloShowerStepMCs, crid, sim, steps);
                  if (diagLevel_ > 1) fillHisto1(cal,sim,info.sims());
              }
              else
              {
                  std::map<SimPtr, std::vector<const StepPointMC*>> newSimStepMap;
                  for (const StepPointMC* step : steps) newSimStepMap[step->simParticle()].push_back(step);
                  for (auto& iter : newSimStepMap)
                  {
                      compressSteps(cal, caloShowerStepMCs, crid, iter.first, iter.second);
                      SimsToKeepUnique.insert(iter.first);
                 }
              }
          }
          ++diagSummary_.nCompressAll_;
          if (doCompress) ++diagSummary_.nCompress_;
      }
      
      //dump the unique set of SimParticles to keep into final vector
      simsToKeep.assign(SimsToKeepUnique.begin(),SimsToKeepUnique.end());

      
      //---------------------------------------------------------------------------------------------------------------
      // Final diag info      
      if (diagLevel_ > 1) 
      {
          hEtot_->Fill(diagSummary_.totalEdep_);
          std::cout<<"CaloShowerStepFromStepPt summary"<<std::endl;
          
          std::set<int> volIds;
          for (auto caloShowerStepMC : caloShowerStepMCs) volIds.insert(caloShowerStepMC.volumeId());
          
          for (auto volId: volIds)
          { 
             std::map<const art::Ptr<SimParticle>, double> simMap;
             for (const auto& caloShowerStepMC : caloShowerStepMCs) if (caloShowerStepMC.volumeId()==volId) simMap[caloShowerStepMC.simParticle()] += caloShowerStepMC.energyDepG4();
             for (auto& kv : simMap) std::cout<<"Vol id: "<<volId<<"  Sim id: "<<kv.first.id()<<"   energy="<<kv.second<<std::endl;
          }
      }      
      
      if (diagLevel_ > 0) 
        std::cout << "[CaloShowerStepFromStepPt::makeCompressedHits] compressed "<<diagSummary_.nCompress_<<" / "<<diagSummary_.nCompressAll_<<" incoming SimParticles"<<std::endl
                  << "[CaloShowerStepFromStepPt::makeCompressedHits] keeping "<<simsToKeep.size()<<" SimParticles"<<std::endl
                  << "[CaloShowerStepFromStepPt::makeCompressedHits] Total sims init: " <<diagSummary_.totalSim_<<std::endl
                  << "[CaloShowerStepFromStepPt::makeCompressedHits] Total caloShower steps: " <<caloShowerStepMCs.size()<<std::endl
                  << "[CaloShowerStepFromStepPt::makeCompressedHits] Total energy deposited / number of stepPointMC: " <<diagSummary_.totalEdep_<<" / "<<diagSummary_.totalStep_<<std::endl
                  << "[CaloShowerStepFromStepPt::makeCompressedHits] Total stepPointMCs seen: " <<diagSummary_.totalChk_<<std::endl;
  }


  //------------------------------------------------------------------------------------------------------------------
  void CaloShowerStepFromStepPt::collectStepBySimAncestor(const Calorimeter& cal, const PhysicalVolumeMultiHelper& vi,
                                                          const HandleVector& stepsHandles, std::map<SimPtr,CaloCompressUtil>& ancestorsMap)
  {
     std::unordered_map<SimPtr,SimPtr> simToAncestorMap;
     for (HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end(); i != e; ++i )
     {
         const art::Handle<StepPointMCCollection>& handle(*i);
         const StepPointMCCollection& steps(*handle);
         for (const auto& step : steps )
         {
             SimPtr sim = step.simParticle();

             SimParticlePtrCollection inspectedSims;
             while (sim->hasParent() && isInsideCalorimeter(cal, vi, sim) )
             {
                 //simparticle starting in one section and ending in another one see note above
                 if (!cal.geomUtil().isContainedSection(sim->startPosition(),sim->endPosition()) ) break;

                 const auto alreadyInspected = simToAncestorMap.find(sim);
                 if (alreadyInspected != simToAncestorMap.end()) {sim = alreadyInspected->second; break;}

                 inspectedSims.push_back(sim);
                 sim = sim->parent();
             }

             for (const SimPtr& inspectedSim : inspectedSims) simToAncestorMap[inspectedSim] = sim;
             ancestorsMap[sim].fill(&step,inspectedSims);

             diagSummary_.totalEdep_ += step.totalEDep();             
         }
         diagSummary_.totalStep_ += steps.size();
      }   
  }





  //-------------------------------------------------------------------------------------------------------------------------
  bool CaloShowerStepFromStepPt::isInsideCalorimeter(const Calorimeter& cal, const PhysicalVolumeMultiHelper& vi, 
                                                     const art::Ptr<SimParticle>& thisSimPtr)
  {
      if (usePhysVol_) return mapPhysVol_.find(&vi.startVolume(*thisSimPtr)) != mapPhysVol_.end();   
      return cal.geomUtil().isInsideCalorimeter(thisSimPtr->startPosition());
  }



  //-----------------------------------------------------------------------------------------------------------------------------------------------
  void CaloShowerStepFromStepPt::collectStepBySim(const HandleVector& stepsHandles, std::map<SimPtr,std::vector<const StepPointMC*>>& simStepMap)
  {
      for ( HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end(); i != e; ++i)
      {
          const art::Handle<StepPointMCCollection>& handle(*i);
          const StepPointMCCollection& steps(*handle);
          for (const auto& step : steps ) simStepMap[step.simParticle()].push_back(&step);
      }
  }


  //-------------------------------------------------------------------------------------------------------------------------------
  void CaloShowerStepFromStepPt::compressSteps(const Calorimeter& cal, CaloShowerStepCollection &caloShowerStepMCs,
                                               int volId, const SimPtr& sim, std::vector<const StepPointMC*>& steps)
  {
     std::sort(steps.begin(), steps.end(), [](const StepPointMC* a, const StepPointMC* b) {return a->time() < b->time();});

     ShowerStepUtil buffer(numZSlices_, ShowerStepUtil::weight_type::energy );

     for (const StepPointMC* step : steps)
     {
         CLHEP::Hep3Vector pos  = cal.geomUtil().mu2eToCrystal(volId,step->position());
         int               idx  = int(std::max(1e-6,pos.z())/zSliceSize_);

         if (buffer.entries(idx)>0 && (step->time()-buffer.t0(idx) > deltaTime_) )
         {
             if (diagLevel_ > 1) {fillHisto2(idx,buffer.energyG4(idx),sim);}
             if (diagLevel_ > 2) {std::cout<<"[CaloShowerStepFromStepPt::compressSteps] inserted  "; buffer.printBucket(idx);}
             diagSummary_.totalChk_ += buffer.entries(idx);
 
             caloShowerStepMCs.push_back(CaloShowerStep(volId, sim, buffer.entries(idx), buffer.time(idx), buffer.energyG4(idx), 
                                                        buffer.energyVis(idx),buffer.pIn(idx),buffer.pos(idx)));
             buffer.reset(idx);
         }

         buffer.add(idx, step->totalEDep(), step->visibleEDep(), step->time(), step->momentum().mag(), pos);
     }

     //do not forget to flush the final buffer(s) :-)
     for (unsigned i=0;i<buffer.nBuckets();++i)
     {
         if (buffer.entries(i) == 0) continue;

         if (diagLevel_ > 1) {fillHisto2(i,buffer.energyG4(i),sim);}
         if (diagLevel_ > 2) {std::cout<<"[CaloShowerStepFromStepPt::compressSteps] inserted ";  buffer.printBucket(i);}
         diagSummary_.totalChk_ += buffer.entries(i);

         caloShowerStepMCs.push_back(CaloShowerStep(volId, sim,  buffer.entries(i), buffer.time(i), buffer.energyG4(i), 
                                                    buffer.energyVis(i),buffer.pIn(i),buffer.pos(i)));
     }
  }

  //-------------------------------------------------------------------------------------------------------------
  bool CaloShowerStepFromStepPt::isCompressible(const int simPdgId, const std::set<int>& processCodes)
  {
      if (compressAll_ || simPdgId > 1000000000)         return true;  //ions are always compressed
      if (procCodes_.find(simPdgId) == procCodes_.end()) return false;

      const std::set<int>& proc = procCodes_[simPdgId];
      for (int code : processCodes) if ( proc.find(code) == proc.end() ) return false;

      return true;
  }

  //-------------------------------------------------------------------------------------------------------------
  void CaloShowerStepFromStepPt::fillHisto1(const Calorimeter& cal, const art::Ptr<SimParticle>& sim,
                                            const std::set<art::Ptr<SimParticle>>& infoSims)
  {
      CLHEP::Hep3Vector startSection = cal.geomUtil().mu2eToDisk(0,sim->startPosition());
      CLHEP::Hep3Vector endSection   = cal.geomUtil().mu2eToDisk(0,sim->endPosition());
      double rStart = sqrt(startSection.x()*startSection.x()+startSection.y()*startSection.y());
      double rEnd   = sqrt(endSection.x()*endSection.x()+endSection.y()*endSection.y());

      hStartPos_->Fill(sim->startPosition().z(),rStart);
      hStopPos_->Fill( sim->endPosition().z(),  rEnd);
      for (const auto& simD: infoSims)
      {
          hStartPos2_->Fill(simD->startPosition().z());
          hStopPos2_->Fill(simD->endPosition().z());
      }
  }
 
  //-------------------------------------------------------------------------------------------------------------
  void CaloShowerStepFromStepPt::fillHisto2(int idx, float edep, const SimPtr& sim)
  {
      hZpos_->Fill(idx);
      hZpos2_->Fill(idx,edep);
      if (sim->genParticle()) hGenId_->Fill(sim->genParticle()->generatorId().id());
  }

  //-------------------------------------------------------------------------------------------------------------
  void CaloShowerStepFromStepPt::dumpAllInfo(const HandleVector& stepsHandles, const Calorimeter& cal)
  {
      std::cout<<"Dumping StepPointMCs  Mu2e / crystal / disk / diskFF frames"<<std::endl;
      for ( HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end(); i != e; ++i )
      {
          const art::Handle<StepPointMCCollection>& handle(*i);
          const StepPointMCCollection& steps(*handle);

          std::cout<<steps.size()<<std::endl;
          for (const auto& step : steps ) 
            std::cout<<step.volumeId()<<" "<<step.totalEDep()<<" "<<step.position()<<" "
                     <<cal.geomUtil().mu2eToCrystal(step.volumeId(),step.position())<<"   "
                     <<cal.geomUtil().mu2eToDisk(cal.crystal(step.volumeId()).diskId(),step.position())<<"   "
                     <<cal.geomUtil().mu2eToDiskFF(cal.crystal(step.volumeId()).diskId(),step.position())<<std::endl;
      }
  }

}

using mu2e::CaloShowerStepFromStepPt;
DEFINE_ART_MODULE(CaloShowerStepFromStepPt);
