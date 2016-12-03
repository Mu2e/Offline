// Transform the deposited energy into scintillator light seen by the crystal readouts
//
// The scintillation light produced by the deposited energy is a random process. It is calculated with a 
// fixed conversion factor fluctuated, and transformed back to energy for convenience. 
//
// The corrections include longitudinal response uniformity, Birke's law (crystal level) and 
// photo-statistics fluctuations (readout level)

// Original author Bertrand Echenard
//

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Selector.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepROCollection.hh"
#include "MCDataProducts/inc/CaloShowerSimCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "SeedService/inc/SeedService.hh"

#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandExponential.h"

#include "TH2D.h"
#include "TH1D.h"

#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <utility>


// Anonymous namespace to hold some helper classes. 
namespace {

     
     class SimParticleEntry
     {
         public:

	     SimParticleEntry(const art::Ptr<mu2e::SimParticle>& sim, const art::Ptr<mu2e::CaloShowerStep>& step, double edepCorr, double time) : 
                sim_(sim),step_(step), edepInit_(step->energyMC()), edepCorr_(edepCorr), time_(time)
             {}

             art::Ptr<mu2e::SimParticle> sim_;
             art::Ptr<mu2e::CaloShowerStep> step_;
             double edepInit_;           
             double edepCorr_;           
             double time_;                
     }; 
    
     
     
     
     class SimParticleSummary
     {
         public:

	     typedef art::Ptr<mu2e::CaloShowerStep> CaloShowerStepPtr;

	     SimParticleSummary(const CaloShowerStepPtr& step, double edepInit, double edepCorr, double time) : 
                steps_(), edepInit_(edepInit), edepCorr_(edepCorr), time_(time), pIn_(step->momentumIn())
             {
		steps_.push_back(step); 
	     }

	     void add(const CaloShowerStepPtr& step, double edepInit, double edepCorr, double time)
	     {
		 steps_.push_back(step);
		 edepInit_ += edepInit;
		 edepCorr_ += edepCorr;
		 time_ = std::min(time_,time);
                 pIn_ = std::max(pIn_,step->momentumIn());
	     }

             std::vector<CaloShowerStepPtr> steps_;
             double edepInit_;           
             double edepCorr_;           
             double time_;
             double pIn_;                
     }; 
} 






namespace mu2e {



  //------------------------------------------------------------
  class CaloShowerStepROFromShowerStep : public art::EDProducer {
 
    public:

        explicit CaloShowerStepROFromShowerStep(fhicl::ParameterSet const& pset) :

          caloShowerStepMCModuleLabel_(pset.get<std::string>("caloShowerStepMCModuleLabel")),
          caloShowerMCName_           (pset.get<std::string>("caloShowerMCName")),
          toff_                       (pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
          blindTime_                  (pset.get<double>     ("blindTime")),
          caloLRUCorrection_          (pset.get<bool>       ("caloLRUcorrection")),
          caloBirksCorrection_        (pset.get<bool>       ("caloBirksCorrection")),
          addTravelTime_              (pset.get<bool>       ("addTravelTime")),
          diagLevel_                  (pset.get<int>        ("diagLevel",0)),
          messageCategory_            ("CaloShowerStepROFromShowerStep"),
          engine_                     ( createEngine(art::ServiceHandle<SeedService>()->getSeed()) ),
          randPoisson_                (engine_),
          randGauss_                  (engine_),          
          randExpo_                   (engine_)          
        {
             produces<CaloShowerStepROCollection>();
             produces<CaloShowerSimCollection>();
        }

        virtual ~CaloShowerStepROFromShowerStep() {}
        virtual void beginJob();
        void produce( art::Event& e);



     private:

         std::string             caloShowerStepMCModuleLabel_;   
         std::string             caloShowerMCName_;   

         SimParticleTimeOffset   toff_; 
         double                  blindTime_;
         double                  mbtime_;

         bool                    caloLRUCorrection_;
         bool                    caloBirksCorrection_;
         bool                    caloPEStatCorrection_;
         bool                    addTravelTime_;

         int                     diagLevel_;
         const std::string       messageCategory_;

         CLHEP::HepRandomEngine& engine_;
         CLHEP::RandPoisson      randPoisson_;
         CLHEP::RandGaussQ       randGauss_;
         CLHEP::RandExponential  randExpo_;

         TH1F*                   hEner_;
         TH2F*                   hLRUCorr_;
         TH1F*                   hBiCorr_;
         TH2F*                   hPECorr_;
         TH1F*                   hPECorr2_;


         void   makeReadoutHits(const art::Handle<CaloShowerStepCollection>&, CaloShowerStepROCollection&, CaloShowerSimCollection&);         
	 double LRUCorrection(int crystalID, double normalizedPosZ, double edepInit, const ConditionsHandle<CalorimeterCalibrations>&);
         double BirksCorrection(int particleCode, double edepInit, const ConditionsHandle<CalorimeterCalibrations>&);
         double photoStatisticsCorrection(int crystalID, double edepInit, double NpePerMeV);
  };


  //----------------------------------------
  void CaloShowerStepROFromShowerStep::beginJob()
  {
      if (diagLevel_ > 2)
      {
         art::ServiceHandle<art::TFileService> tfs;
         hEner_      = tfs->make<TH1F>("hEner",   "energy deposited / hit",   200,    0,  20);
         hLRUCorr_   = tfs->make<TH2F>("hLRUCorr","LRU correction",           100,    0,   1, 100, -0.5, 0.5);
         hBiCorr_    = tfs->make<TH1F>("hBiCorr", "Birkes stat correction",   100,    0,   1);
         hPECorr_    = tfs->make<TH2F>("hPECorr", "PE stat correction",       100,    0., 50, 100,    0,  50);
         hPECorr2_   = tfs->make<TH1F>("hPECorr2","PE stat correction",       100,  -0.5,    0.5);
     }
  }





  void  CaloShowerStepROFromShowerStep::produce(art::Event& event) 
  {
      if ( diagLevel_ > 0) std::cout << "[CaloShowerStepROFromShowerStep::produce] begin" << std::endl;
      
      art::ServiceHandle<GeometryService> geom;    
      if (!geom->hasElement<Calorimeter>()) return;


      //update condition cache
      ConditionsHandle<AcceleratorParams> accPar("ignored");
      mbtime_ = accPar->deBuncherPeriod;
      toff_.updateMap(event);


      // A container to hold the output hits.
      std::unique_ptr<CaloShowerStepROCollection> caloShowerStepROs(new CaloShowerStepROCollection);
      
      // A container to hold the output hits.
      std::unique_ptr<CaloShowerSimCollection> caloShowerSims(new CaloShowerSimCollection);


      //Get handles to caloStepMC collection
      art::Handle<CaloShowerStepCollection> caloShowerStepMCHandle;
      event.getByLabel(caloShowerStepMCModuleLabel_, caloShowerMCName_,   caloShowerStepMCHandle);


      makeReadoutHits(caloShowerStepMCHandle, *caloShowerStepROs, *caloShowerSims);

      // Add the output hit collection to the event
      event.put(std::move(caloShowerStepROs));
      event.put(std::move(caloShowerSims));
      
      if ( diagLevel_ > 0) std::cout << "[CaloShowerStepROFromShowerStep::produce] end" << std::endl;
 } 








//------------------------------------------------------------------------------------------------------------------------------------- 
  void CaloShowerStepROFromShowerStep::makeReadoutHits(const art::Handle<CaloShowerStepCollection>& caloShowerStepMCCollHandle, 
                                                       CaloShowerStepROCollection& caloShowerStepROs, CaloShowerSimCollection& caloShowerSims)
  {

       GlobalConstantsHandle<ParticleDataTable>  pdt;
       ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

       const Calorimeter& cal  = *(GeomHandle<Calorimeter>());
       int    nROs             = cal.caloGeomInfo().nROPerCrystal();
       double cryhalflength    = cal.caloGeomInfo().crystalHalfLength();
       double refractiveIndex  = cal.caloGeomInfo().refractiveIndex();
       double lightSpeed       = 300; // mm/ns


       const CaloShowerStepCollection& caloShowerSteps(*caloShowerStepMCCollHandle);
       const CaloShowerStep*           caloShowerStepMCBase = &caloShowerSteps.front();

       std::map<int, std::vector<SimParticleEntry> > simEntriesMap;


       //-----------------------------------------------------------------------
       //store corrected energy deposits for each redouts

       double totalEdep(0.0),totalEdepCorr(0.0);
       int    totalSteps(0),totalHit(0);

       for (auto istep = caloShowerSteps.begin(); istep !=caloShowerSteps.end(); ++istep)
       {        
           const CaloShowerStep& step = *istep;

	   // time folding and filtering, see docdb-3425 for a stunning explanation
           double hitTimeUnfolded = toff_.totalTimeOffset(istep->simParticle())+istep->timeStepMC();
           double hitTime         = fmod(hitTimeUnfolded,mbtime_);

           if (hitTime < blindTime_ || hitTime > mbtime_ ) continue;

           size_t idx = (&step - caloShowerStepMCBase);
           art::Ptr<CaloShowerStep> stepPtr = art::Ptr<CaloShowerStep>(caloShowerStepMCCollHandle,idx);

           int    crystalID = step.volumeId();
	   int    ROIDBase  = cal.ROBaseByCrystal(crystalID);
           int    pdgId     = step.simParticle()->pdgId();
           double edep_corr = step.energyMC();
           double posZ      = step.position().z();
	   double timeToRO  = (2.0*cryhalflength-posZ)*refractiveIndex/lightSpeed;

	   const  art::Ptr<SimParticle> sim = step.simParticle();


	   // apply travel time to readout correction 
	   if (addTravelTime_) hitTime += timeToRO;


           // apply corrections on energy deposition -> scintillation light in the crystal, see Note
           if (caloLRUCorrection_)
                edep_corr = LRUCorrection(crystalID, posZ/cryhalflength, edep_corr, calorimeterCalibrations);

           if (caloBirksCorrection_)
                edep_corr = BirksCorrection(pdgId,edep_corr, calorimeterCalibrations);


	   // energy equivalent deposited in each readout 	    
	   for (int i=0; i<nROs; ++i)
	   {
               int ROID = ROIDBase + i;

               //apply photo-statistics correction 
               double NpePerMeV = calorimeterCalibrations->peMeV(ROID);

               double edep_corr_RO(edep_corr);
	       if (caloPEStatCorrection_) edep_corr_RO = photoStatisticsCorrection(ROID, edep_corr_RO, NpePerMeV);          
               if (edep_corr_RO < 1e-4) continue;

	       //finally, add a caloShowerStepRO entry here 
	       caloShowerStepROs.push_back(CaloShowerStepRO(ROID,stepPtr,hitTime,edep_corr_RO));
	   }

	   simEntriesMap[crystalID].push_back(SimParticleEntry(sim,stepPtr,edep_corr,hitTime));


           if (diagLevel_ > 0) 
	   {
	        totalEdep += step.energyMC(); 
		totalEdepCorr += edep_corr; 
		totalSteps += step.nCompress();
		++totalHit;          

		if (diagLevel_ > 2) hEner_->Fill(edep_corr); 
	        if (diagLevel_ > 3) std::cout<<"[CaloShowerStepROFromShowerStep::produce] "<<crystalID<<" "<<step.simParticle()<<" "<<stepPtr<<" "<<hitTime<<" "<<edep_corr<<std::endl;
                if (diagLevel_ > 3) std::cout<<"   "<<step.simParticle()->startPosition()<<" "<<step.simParticle()->endPosition()<<" "<<  step.simParticle()->pdgId()<<std::endl; 	    
           }	            
       }


       
       //-----------------------------------------------------------------------
       // Produce the CaloShowerSim from the simEntriesMap
       
       for (auto& kv : simEntriesMap)
       {	   
	   int crystalID = kv.first;
	   std::vector<SimParticleEntry> entries = kv.second;

	   // fill the summary map for each simPtr for a given crystalID
	   std::map< art::Ptr<SimParticle>, SimParticleSummary> simSumMap;
	   for (auto& entry : entries)
	   {
	       auto mfind = simSumMap.find(entry.sim_);
	       if (mfind==simSumMap.end()) 
		   simSumMap.insert( std::make_pair(entry.sim_,SimParticleSummary(entry.step_,entry.edepInit_,entry.edepCorr_,entry.time_)) ); 
	       else 
		   mfind->second.add(entry.step_,entry.edepInit_,entry.edepCorr_,entry.time_);   
	   }

	   // create the CaloShowerSim entries	   
	   for (auto& kv2 : simSumMap)
	   {
	       auto simPtr  = kv2.first;
	       auto summary = kv2.second;
	       caloShowerSims.push_back(CaloShowerSim(crystalID,simPtr,summary.steps_, summary.time_,summary.edepCorr_,summary.edepInit_, summary.pIn_));
	   }
	   
       }
         


       //-----------------------------------------------------------------------
       // diagnosis

       if (diagLevel_ > 0) std::cout<<"[CaloShowerStepROFromShowerStep::produce] found energy (energy corr) / nStepsMC / nCORHit "
                                    <<totalEdep<<" ("<<totalEdepCorr<<") / "<<totalSteps<<" / "<<totalHit<<std::endl;

       if (diagLevel_ > 3)
       {
            std::cout<<"Checking Sims"<<std::endl;
            double csmEtot(0);
            for (auto& csm :  caloShowerSims) 
            {          
                csmEtot +=csm.energyMC();
                std::cout<<csm.crystalId()<<" "<<csm.sim()<<" "<<csm.time()<<" "<<csm.energy()<<" "<<csm.energyMC()<<std::endl;
                for (auto& st : csm.caloShowerSteps()) std::cout<<"  "<<st<<std::endl;
            }

            std::cout<<"CSM Etot "<<csmEtot<<std::endl;
       }

  }





  

 

  //----------------------------------------------------------------------------------------------------------------------------------
  // apply a correction of type Energy = ((1-s)*Z/HL+s)*energy where Z position along the crystal, HL is the crystal half-length
  // and s is the intercept at Z=0 (i.e. non-uniformity factor, e.g. 5% -> s = 1.05)
  
  double CaloShowerStepROFromShowerStep::LRUCorrection(int crystalID, double normalizedPosZ, double edepInit, const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
      double edep(edepInit);     
     
      double alpha  = calorimeterCalibrations->LRUpar0(crystalID);
      double factor =  (1.0-alpha)*normalizedPosZ +alpha;
      edep *= factor;

      if (diagLevel_ > 2) hLRUCorr_->Fill(normalizedPosZ,edep/edepInit-1);
      if (diagLevel_ > 3) std::cout<<"[CaloShowerStepROFromShowerStep::LRUCorrection] before / after LRU -> edep_corr = "<< edep<<"  /  "<<edepInit<<"  at position Z="<<normalizedPosZ<<std::endl;          
      return edep;     
  }

  //-----------------------------------------------------------------------------
  double CaloShowerStepROFromShowerStep::BirksCorrection(int particleCode, double edepInit,const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
      double edep(edepInit);
      
      if (particleCode==2212 || particleCode == 2112) edep /= calorimeterCalibrations->BirkCorrHadron();      
      
      if (diagLevel_ > 2) hBiCorr_->Fill(edep/edepInit);
      return edep;    
  }

  //-------------------------------------------------------------------------------------------------------------------------------------------------------
  double CaloShowerStepROFromShowerStep::photoStatisticsCorrection(int crystalID, double edepInit, double NpePerMeV)
  { 
      double lightYield = randPoisson_.fire(edepInit*NpePerMeV);
      double edep       = lightYield/NpePerMeV;
      
      if (diagLevel_ > 2) hPECorr_->Fill(edep,edepInit);
      if (diagLevel_ > 2 && lightYield>10) hPECorr2_->Fill((edep-edepInit)/sqrt(edepInit));            
      return edep;
  }


} 

using mu2e::CaloShowerStepROFromShowerStep;
DEFINE_ART_MODULE(CaloShowerStepROFromShowerStep);




/*
//little snipper I kept for convenience

double crystalDecayTime = cal.caloGeomInfo().crystalDecayTime();

////generate each photon hitting the readout, correcting time for scintillator decay time - get crackin'
//int nPhot            = std::max(int(edep_corr_RO*NpePerMeV),1);
//double energyPerPhot = edep_corr_RO/float(nPhot);

//std::vector<double> photonTime(nPhot);
//std::generate(photonTime.begin(),photonTime.end(),[&](){return hitTime+randExpo_.fire(crystalDecayTime);});  
*/
