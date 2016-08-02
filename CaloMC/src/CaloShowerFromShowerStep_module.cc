//
// An EDProducer Module that reads CaloShowerStepPointMC and produces CaloHit, CaloHitMCTruth objects
//
// The scintillation light produced by the deposited energy is a random process, and depends on the position and other factors.
// The scintillation light is calculated with a fixed conversion factor, then fluctuated, and transformed back to energy for convenience. 
// The corrections include longitudinal response uniformity, non-linearity and Birke's law.
// The photo-statistics and electronic noise are simulated for each sensor independently.

// Original author Bertrand Echenard
//
// C++ includes.
#include <iostream>
#include <string>
#include <cmath>
#include <map>
#include <vector>
#include <utility>

// Framework includes.
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

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerCollection.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "SeedService/inc/SeedService.hh"


// Other includes.
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "TH2D.h"
#include "TH1D.h"




namespace mu2e {



    // Anonymous namespace to hold some helper classes. COR means Crystal Or Readout
    namespace {

         class CORHit {

             public:
               
               CORHit(art::Ptr<CaloShowerStep> const& step, double edep_corr,double edep_uncorr, double time) : 
                 step_(step), edep_corr_(edep_corr), edep_uncorr_(edep_uncorr),time_(time)
               {}

               art::Ptr<CaloShowerStep> step_;
               double                   edep_corr_;           
               double                   edep_uncorr_;           
               double                   time_;                
         }; 
    } 








  //------------------------------------------------------------
  class CaloShowerFromShowerStep : public art::EDProducer {
 
    public:

        explicit CaloShowerFromShowerStep(fhicl::ParameterSet const& pset) :

          caloShowerStepMCModuleLabel_(pset.get<std::string>("caloShowerStepMCModuleLabel")),
          caloShowerMCName_           (pset.get<std::string>("caloShowerMCName")),
          caloLRUCorrection_          (pset.get<bool>       ("caloLRUcorrection")),
          caloBirkesCorrection_       (pset.get<bool>       ("caloBirkesCorrection")),
          toff_                       (pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
          blindTime_                  (pset.get<double>     ("blindTime")),
          timeGap_                    (pset.get<double>     ("timeGap")),
          filterEnergy_               (pset.get<double>     ("filterEnergy")),
          filterDeltaT_               (pset.get<double>     ("filterDeltaT")),
          diagLevel_                  (pset.get<int>        ("diagLevel",0)),
          messageCategory_            ("CaloShowerFromShowerStep"),
          engine_                     ( createEngine(art::ServiceHandle<SeedService>()->getSeed()) ),
          randGauss_                  (engine_)          
        {
             produces<CaloShowerCollection>();
        }

        virtual ~CaloShowerFromShowerStep() { }
        virtual void beginJob();
        void produce( art::Event& e);



     private:

         std::string             caloShowerStepMCModuleLabel_;   
         std::string             caloShowerMCName_;   

         bool                    caloPEStatCorrection_;
         bool                    caloLRUCorrection_;
         bool                    caloBirkesCorrection_;

         SimParticleTimeOffset   toff_; 
         double                  blindTime_;
         double                  mbtime_;
         double                  timeGap_;
         
         double                  filterEnergy_;
         double                  filterDeltaT_;

         int                     diagLevel_;
         const std::string       messageCategory_;

         CLHEP::HepRandomEngine& engine_;
         CLHEP::RandGaussQ       randGauss_;

         TH1F*                   hEner_;
         TH1F*                   hinit_;
         TH1F*                   hfinal_;
         TH2F*                   hLRUCorr_;
         TH1F*                   hBiCorr_;

         void   makeCalorimeterHits (const art::Handle<CaloShowerStepCollection>&, CaloShowerCollection &);
         double LRUCorrection(double normalizedPosZ, double energy,int crystalId, const ConditionsHandle<CalorimeterCalibrations>&);
         double BirkesCorrection(int pdgId, double energy);



  };

  //----------------------------------------
  void CaloShowerFromShowerStep::beginJob()
  {
      art::ServiceHandle<art::TFileService> tfs;
      hEner_      = tfs->make<TH1F>("hEner",   "energy deposited / hit",        200,    0,  20);
      hinit_      = tfs->make<TH1F>("hInit",   "Initial number of calo showers", 25,    0,  25);
      hfinal_     = tfs->make<TH1F>("hFinal",  "Final number of calo showers",   25,    0,  25);
      hLRUCorr_   = tfs->make<TH2F>("hLRUCorr","LRU correction",                100,    0, 200, 100, -0.5, 0.5);
      hBiCorr_    = tfs->make<TH1F>("hBiCorr", "Birkes stat correction",        100,    0,   1);
 }





  void  CaloShowerFromShowerStep::produce(art::Event& event) 
  {
      art::ServiceHandle<GeometryService> geom;    
      if( !(geom->hasElement<Calorimeter>()) ) return;


      //update condition cache
      ConditionsHandle<AcceleratorParams> accPar("ignored");
      mbtime_ = accPar->deBuncherPeriod;
      toff_.updateMap(event);


      // A container to hold the output hits.
      std::unique_ptr<CaloShowerCollection>  caloShowers(new CaloShowerCollection);


      //Get handles to caloStepMC collection
      art::Handle<CaloShowerStepCollection> caloShowerStepMCHandle,caloROShowerStepMCHandle;
      event.getByLabel(caloShowerStepMCModuleLabel_, caloShowerMCName_,   caloShowerStepMCHandle);


      makeCalorimeterHits(caloShowerStepMCHandle, *caloShowers);

      // Add the output hit collection to the event
      event.put(std::move(caloShowers));
  } 








//------------------------------------------------------------------------------------------------------------------------------------- 
  void CaloShowerFromShowerStep::makeCalorimeterHits(const art::Handle<CaloShowerStepCollection>& caloShowerStepMCCollHandle, 
                                                     CaloShowerCollection &caloShowers)
  {

       GlobalConstantsHandle<ParticleDataTable>  pdt;
       ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

       const Calorimeter& cal = *(GeomHandle<Calorimeter>());
       double cryhalflength   = cal.caloGeomInfo().crystalHalfLength();


       const CaloShowerStepCollection& caloShowerSteps(*caloShowerStepMCCollHandle);
       const CaloShowerStep*           caloShowerStepMCBase = &caloShowerSteps.front();



       //-----------------------------------------------------------------------
       //store corrected energy deposits for each crystal hits

       double totalEdep(0.0),totalEdepCorr(0.0),eDepCorrTot(0);
       int    totalSteps(0), totalCORHit(0), numpreFilter(0), numpostFilter(0),nCorrTot(0);

       std::map<int, std::vector<CORHit> >  hitMapCrystal;    

       for (auto istep = caloShowerSteps.begin(); istep !=caloShowerSteps.end(); ++istep)
       {        
            const CaloShowerStep& step = *istep;

            size_t idx = (&step - caloShowerStepMCBase);
            art::Ptr<CaloShowerStep> stepPtr = art::Ptr<CaloShowerStep>(caloShowerStepMCCollHandle,idx);

            int crystalId    = step.volumeId();
            int pdgId        = step.simParticle()->pdgId();
            double edep_corr = step.energyMC();
            double posZ      = step.position().z();

            if (diagLevel_ > 0) {totalEdep += step.energyMC(); totalSteps += step.nCompress();}

            //apply corrections on energy deposition -> scintillation light, see Note

            if (caloLRUCorrection_)
                 edep_corr = LRUCorrection(posZ/cryhalflength, edep_corr, crystalId, calorimeterCalibrations);

            if (caloBirkesCorrection_)
                 edep_corr = BirkesCorrection(pdgId,edep_corr);


            if (diagLevel_ > 2)
            {
                double edep_init = step.energyMC();
                hLRUCorr_->Fill(posZ,LRUCorrection(posZ/cryhalflength, edep_init, crystalId, calorimeterCalibrations)/edep_init-1);
                hBiCorr_->Fill(BirkesCorrection(pdgId,edep_init)/edep_init);    
            } 


            // time folding, see docdb-3425 for a stunning explanation
            double hitTimeUnfolded = toff_.totalTimeOffset(step.simParticle())+step.timeStepMC();
            double hitTime         = fmod(hitTimeUnfolded,mbtime_);

            if (hitTime < blindTime_) continue;

            hitMapCrystal[crystalId].emplace_back(CORHit(stepPtr,edep_corr,step.energyMC(),hitTime));
            
            if (diagLevel_ > 0) {totalEdepCorr += edep_corr; ++totalCORHit;}          
       }



       if (diagLevel_ > 0) std::cout<<"CaloShowerFromShowerStep found energy (energy corr) / nStepsMC / nCORHit "
                                   <<totalEdep<<" ("<<totalEdepCorr<<") / "<<totalSteps<<" / "<<totalCORHit<<std::endl;




       //Group hits within short time window
       for (auto &kv : hitMapCrystal )
       {
           int crystalId = kv.first;
           std::vector<CORHit> &hits = kv.second;
           std::vector<CaloShower> caloShowerBuild;

           std::sort(hits.begin(),hits.end(), [](const CORHit& a, const CORHit& b){return a.time_ < b.time_;} );

           double h_time        = hits[0].time_;
           double h_edep        = hits[0].edep_corr_;
           double h_edep_uncorr = hits[0].edep_uncorr_;
           std::vector<art::Ptr<CaloShowerStep>> h_sims {hits[0].step_};

           for( size_t i=1; i<hits.size(); ++i )
           {          
               //Save hit and create new hits
               if ( (hits[i].time_ - h_time) > timeGap_ )
               {            
                   if (diagLevel_ > 2) hEner_->Fill(h_edep);
                   if (diagLevel_ > 0) {eDepCorrTot+= h_edep; nCorrTot += h_sims.size(); }                   

                   caloShowerBuild.push_back(CaloShower(crystalId,h_sims,h_time,h_edep,h_edep_uncorr));
                   h_sims.clear();
                   h_time        = hits[i].time_;
                   h_edep        = hits[i].edep_corr_;
                   h_edep_uncorr = hits[i].edep_uncorr_;
                   h_sims.push_back(hits[i].step_);
               }           
               else
               {
                   h_edep += hits[i].edep_corr_;
                   h_edep_uncorr += hits[i].edep_uncorr_;
                   h_sims.push_back(hits[i].step_);
               }
           }

           //do not forget to save the last hit!
           caloShowerBuild.push_back(CaloShower(crystalId,h_sims,h_time,h_edep,h_edep_uncorr));


           if (diagLevel_ > 0)
           {
               eDepCorrTot += h_edep;
               nCorrTot += h_sims.size();
               hinit_->Fill(caloShowerBuild.size());
               numpreFilter += caloShowerBuild.size();
           }

           
           //now a little bit of pre-filtering, remove low-energy isolated hits
           if (filterEnergy_ > 1e-3)
           {
               auto iterBefore  = caloShowerBuild.rend();
               auto iterCurrent = caloShowerBuild.begin();
               auto iterAfter   = caloShowerBuild.begin();

               while(iterCurrent != caloShowerBuild.end())
               {
                   ++iterAfter;     

                   double deltaTimePlus(1e6), deltaTimeMinus(1e6);
                   if (iterAfter   != caloShowerBuild.end() )  deltaTimePlus  = iterAfter->time() - iterCurrent->time();
                   if (iterBefore  != caloShowerBuild.rend())  deltaTimeMinus = iterCurrent->time() - iterBefore->time();

                   if (iterCurrent->energy() < filterEnergy_ && std::min(deltaTimePlus,deltaTimeMinus) > filterDeltaT_)
                   {
                       iterCurrent = caloShowerBuild.erase(iterCurrent);                      
                   } 
                   else
                   {
                       ++iterBefore;
                       ++iterCurrent;                
                   } 
               }
           }
           //finally copy the remaining caloshowers to the final vector 
           std::copy(caloShowerBuild.begin(),caloShowerBuild.end(),std::back_inserter(caloShowers));

           
           if (diagLevel_ > 0) 
           {              
              numpostFilter += caloShowerBuild.size();
              if (diagLevel_ > 2) hfinal_->Fill(caloShowerBuild.size());
           }

      }


      if (diagLevel_ > 0)
      {
         std::cout<<"CaloShowerFromShowerStep produced caloShower prefilter / postfilter "<<numpreFilter<<" / "<<numpostFilter<<std::endl;
         std::cout<<"CaloShowerFromShowerStep corrected energy / nCOR tot (prefiltre)  "<<eDepCorrTot<<" / "<<nCorrTot<<std::endl;
      } 


  }





  

 

  //----------------------------------------------------------------------------------------------------------------------------------
  // apply a correction of type Energy = ((1-s)*Z/HL+s)*energy where Z position along the crystal, HL is the crystal half-length
  // and s is the intercept at Z=0 (i.e. non-uniformity factor, e.g. 5% -> s = 1.05)
  
  double CaloShowerFromShowerStep::LRUCorrection(double normalizedPosZ, double energy, int crystalId, const ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
      double edep(energy);     
     
      double alpha  = calorimeterCalibrations->LRUpar0(crystalId);
      double factor =  (1.0-alpha)*normalizedPosZ +alpha;
      edep *= factor;

      if (diagLevel_ > 3) std::cout<<"CaloShowerFromShowerStep before / after LRU -> edep_corr = "<< edep<<"  /  "<<energy<<"  at position Z="<<normalizedPosZ<<std::endl;          
      return edep;     
  }



  //---------------------------------------------------------------------------------------
  double CaloShowerFromShowerStep::BirkesCorrection(int particleCode, double energy)
  {
      double edep(energy);
      if (particleCode==2212 || particleCode == 2112) edep /= 4.0;
      return edep;    
  }




} 

using mu2e::CaloShowerFromShowerStep;
DEFINE_ART_MODULE(CaloShowerFromShowerStep);


