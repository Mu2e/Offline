//
// An EDProducer Module that reads StepPointMC and produces CaloShower objects
//
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
#include "MCDataProducts/inc/CaloShowerCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
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








  //--------------------------------------------------------------------
  class CaloShowerFromStepPt : public art::EDProducer {
 
    public:

        explicit CaloShowerFromStepPt(fhicl::ParameterSet const& pset) :

          // Parameters
          caloStepPtsModuleLabel_ (pset.get<std::string>("caloStepPtsModuleLabel")),
          caloLRUCorrection_      (pset.get<bool>       ("caloLRUcorrection")),
          caloBirkesCorrection_   (pset.get<bool>       ("caloBirkesCorrection")),
          toff_                   (pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
          blindTime_              (pset.get<double>     ("blindTime")),
          timeGap_                (pset.get<double>     ("timeGap")),
          filterEnergy_           (pset.get<double>     ("filterEnergy")),
          filterDeltaT_           (pset.get<double>     ("filterDeltaT")),
          diagLevel_              (pset.get<int>        ("diagLevel",0)),
          messageCategory_        ("CaloShowerFromStepPt"),
          engine_                 ( createEngine(art::ServiceHandle<SeedService>()->getSeed()) ),
          randGauss_              (engine_)          
        {
             produces<CaloShowerCollection>();
        }

        virtual ~CaloShowerFromStepPt() { }
        virtual void beginJob();
        void produce( art::Event& e);



     private:

         typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;
     
         std::string             caloStepPtsModuleLabel_;   

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
         TH2F*                   hLRUCorr_;
         TH1F*                   hBiCorr_;

         void   makeCalorimeterHits(const HandleVector& crystalStepsHandles, CaloShowerCollection &caloShowers);
         double LRUCorrection(double normalizedPosZ, double energy,int crystalId, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);
         double BirkesCorrection(int pdgId, double energy, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);



  };


  //-----------------------------------------------
  void CaloShowerFromStepPt::beginJob()
  {
      if (diagLevel_ > 2)
      {
         art::ServiceHandle<art::TFileService> tfs;
         hEner_      = tfs->make<TH1F>("hEner",   "energy deposited / hit",     200,    0,  20);
         hLRUCorr_   = tfs->make<TH2F>("hLRUCorr","LRU correction",             100,    0, 200, 100, -0.5, 0.5);
         hBiCorr_    = tfs->make<TH1F>("hBiCorr", "Birkes stat correction",     100,    0,   1);
      }
  }





  //------------------------------------------------------------
  void CaloShowerFromStepPt::produce(art::Event& event) 
  {
      art::ServiceHandle<GeometryService> geom;    
      if( !(geom->hasElement<Calorimeter>()) ) return;


      //update condition cache
      ConditionsHandle<AcceleratorParams> accPar("ignored");
      mbtime_ = accPar->deBuncherPeriod;
      toff_.updateMap(event);


      // These selectors will select data products with the given instance name, and ignore all other fields of the product ID.
      art::ProductInstanceNameSelector getCrystalSteps(caloStepPtsModuleLabel_);
      HandleVector crystalStepsHandles;      
      event.getMany( getCrystalSteps, crystalStepsHandles);


      std::unique_ptr<CaloShowerCollection>  caloShowers(new CaloShowerCollection);

      makeCalorimeterHits(crystalStepsHandles, *caloShowers);

      event.put(std::move(caloShowers));
  } 








  //------------------------------------------------------------------------------------------------------------------------- 
  void CaloShowerFromStepPt::makeCalorimeterHits(const HandleVector& crystalStepsHandles, CaloShowerCollection &caloShowers)
  {

       GlobalConstantsHandle<ParticleDataTable>  pdt;
       ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

       const Calorimeter& cal = *(GeomHandle<Calorimeter>());
       double cryhalflength   = cal.caloGeomInfo().crystalHalfLength();


       //-----------------------------------------------------------------------
       //store corrected energy deposits for each crystal hits

       double totalEdep(0.0),totalEdepCorr(0.0),eDepCorrTot(0);
       int    totalSteps(0), totalCORHit(0),nCorrTot(0);

       std::map<int, std::vector<CORHit> >  hitMapCrystal;           
       art::Ptr<CaloShowerStep> stepPtr;
       
       
       
       for ( HandleVector::const_iterator i=crystalStepsHandles.begin(), e=crystalStepsHandles.end(); i != e; ++i )
       {     
	    const art::Handle<StepPointMCCollection>& handle(*i);
	    const StepPointMCCollection& stepMCs(*handle);

	    for (auto istep = stepMCs.begin(); istep !=stepMCs.end(); ++istep)
	    {        
		  const StepPointMC& step(*istep);
	          CLHEP::Hep3Vector pos  = cal.toCrystalFrame(step.volumeId(),step.position());

        	  int crystalId    = step.volumeId();
        	  int pdgId        = step.simParticle()->pdgId();
        	  double edep_corr = step.eDep();		  
        	  double posZ      = pos.z();

        	  if (diagLevel_ > 0) {totalEdep += edep_corr; totalSteps += 1;}

        	  //apply corrections on energy deposition -> scintillation light, see Note
        	  if (caloLRUCorrection_)
                       edep_corr = LRUCorrection(posZ/cryhalflength, edep_corr, crystalId, calorimeterCalibrations);

        	  if (caloBirkesCorrection_)
                       edep_corr = BirkesCorrection(pdgId,edep_corr, calorimeterCalibrations);


        	  if (diagLevel_ > 2)
        	  {
                      double edep_init = step.eDep();
                      hLRUCorr_->Fill(posZ,LRUCorrection(posZ/cryhalflength, edep_init, crystalId, calorimeterCalibrations)/edep_init-1);
                      hBiCorr_->Fill(BirkesCorrection(pdgId,edep_init,calorimeterCalibrations)/edep_init);    
        	  } 


        	  // time folding, see docdb-3425 for a stunning explanation
        	  double hitTimeUnfolded = toff_.totalTimeOffset(step.simParticle())+step.time();
        	  double hitTime         = fmod(hitTimeUnfolded,mbtime_);

        	  if (hitTime < blindTime_) continue;

        	  hitMapCrystal[crystalId].emplace_back(CORHit(stepPtr,edep_corr,step.eDep(),hitTime));

        	  if (diagLevel_ > 0) {totalEdepCorr += edep_corr; ++totalCORHit;}          
	     }
	}




       if (diagLevel_ > 0) std::cout<<"CaloShowerFromStepPt found energy (energy corr) / nStepsMC / nCORHit "
                                   <<totalEdep<<" ("<<totalEdepCorr<<") / "<<totalSteps<<" / "<<totalCORHit<<std::endl;




       //Group hits within short time window
       for (auto &kv : hitMapCrystal )
       {
           int crystalId = kv.first;
           std::vector<CORHit> &hits = kv.second;

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

                   caloShowers.push_back(CaloShower(crystalId,h_sims,h_time,h_edep,h_edep_uncorr));
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
           caloShowers.push_back(CaloShower(crystalId,h_sims,h_time,h_edep,h_edep_uncorr));

      }


      if (diagLevel_ > 0)
      {
         std::cout<<"CaloShowerFromStepPt corrected energy / nCOR tot (prefiltre)  "<<eDepCorrTot<<" / "<<nCorrTot<<std::endl;
         std::cout<<"CaloShowerFromStepPt caloShower size "<< caloShowers.size()<<std::endl;
      } 


  }





  

 

  //-----------------------------------------------------------------------------------------------------------------------------
  // apply a correction of type Energy = ((1-s)*Z/HL+s)*energy where Z position along the crystal, HL is the crystal half-length
  // and s is the intercept at Z=0 (i.e. non-uniformity factor, e.g. 5% -> s = 1.05)
  
  double CaloShowerFromStepPt::LRUCorrection(double normalizedPosZ, double energy, int crystalId, 
                                             ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
      double edep(energy);     
     
      double alpha  = calorimeterCalibrations->LRUpar0(crystalId);
      double factor =  (1.0-alpha)*normalizedPosZ +alpha;
      edep *= factor;

      if (diagLevel_ > 3) std::cout<<"CaloShowerFromStepPt before / after LRU -> edep_corr = "<< edep<<"  /  "<<energy<<"  at position Z="<<normalizedPosZ<<std::endl;          
      return edep;     
  }



  //----------------------------------------------------------------------------------------------------------------
  double CaloShowerFromStepPt::BirkesCorrection(int particleCode, double energy, 
                                                ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
      double edep(energy);      
      if (particleCode==2212 || particleCode == 2112) edep /= calorimeterCalibrations->BirkCorrHadron();      
      return edep;    
  }




} 

using mu2e::CaloShowerFromStepPt;
DEFINE_ART_MODULE(CaloShowerFromStepPt);


