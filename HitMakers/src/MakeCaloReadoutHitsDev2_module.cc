// To complete: corrections!
//
// An EDProducer Module that reads CaloShowerStepPointMC and produces CaloHit, CaloHitMCTruth objects
//
// Original author Bertrand Echenard
//
// Note: there is still a lot of work to do on the corrections...

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
#include "CalorimeterGeom/inc/sort_functors.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "DataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/CaloShowerStepMCCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepMC.hh"
#include "MCDataProducts/inc/CaloShowerCollection.hh"
#include "MCDataProducts/inc/CaloShower.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "SeedService/inc/SeedService.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "TH1D.h"




namespace mu2e {



    // Anonymous namespace to hold some helper classes. COR is an acronym for Crystal Or Readout
    namespace {

	 class CORHit {

	     public:
               
	       CORHit(art::Ptr<CaloShowerStepMC> const& step, double edep_corr, double time) : 
	         step_(step), edep_corr_(edep_corr), time_(time)
	       {}

	       art::Ptr<CaloShowerStepMC> step_;
	       double                     edep_corr_;           
	       double                     time_;                
	 }; 
    } 








  //--------------------------------------------------------------------
  class MakeCaloReadoutHitsDev2 : public art::EDProducer {
 
    public:

	explicit MakeCaloReadoutHitsDev2(fhicl::ParameterSet const& pset) :

	  // Parameters
	  _caloShowerStepMCModuleLabel(pset.get<std::string>("caloShowerStepMCModuleLabel")),
	  _caloShowerMCName(           pset.get<std::string>("caloShowerMCName")),
	  _caloLRUcorrection(          pset.get<bool>("caloLRUcorrection",0)),
	  _caloNonLinCorrection(       pset.get<bool>("caloNonLinCorrection",0)),
	  _caloBirkesCorrection(       pset.get<bool>("caloBirkesCorrection",0)),
	  _toff(                       pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
	  _blindTime(                  pset.get<double>("blindTime",300.)),
	  _timeGap(                    pset.get<double>("timeGap",150.)),
	  _filterEnergy(               pset.get<double>("filterEnergy",0.3)),
	  _filterDeltaT(               pset.get<double>("filterDeltaT",150)),
	  _messageCategory(            "CaloReadoutHitsMakerDev2"),
	  _diagLevel(                  pset.get<int>("diagLevel",0)),
	  _maxFullPrint(               pset.get<int>("maxFullPrint",5)),
	  _randGauss(                  createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
	  _hEtot(0),_hEner(0),_hinit(0),_hfinal(0)
	{
	  produces<CaloShowerCollection>();
	}

	virtual ~MakeCaloReadoutHitsDev2() { }
	virtual void beginJob();
	void produce( art::Event& e);



     private:

	 std::string            _caloShowerStepMCModuleLabel;   
	 std::string            _caloShowerMCName;   

	 bool                   _caloLRUcorrection;
	 bool                   _caloNonLinCorrection;
	 bool                   _caloBirkesCorrection;

	 SimParticleTimeOffset  _toff; 
	 double                 _blindTime;
	 double                 _mbtime;
	 double                 _timeGap;
	 
	 double                 _filterEnergy;
	 double                 _filterDeltaT;

	 const std::string      _messageCategory;
	 int                    _diagLevel;  
	 int                    _maxFullPrint;
	 CLHEP::RandGaussQ      _randGauss;

	 TH1F*                  _hEtot;
	 TH1F*                  _hEner;
	 TH1F*                  _hinit;
	 TH1F*                  _hfinal;




	 void makeCalorimeterHits (art::Handle<CaloShowerStepMCCollection>const&, 
                        	   CaloShowerCollection &);

	 double nonLinearityCorrection(int cryId, int particleCode, double energy, 
				       ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations,
				       GlobalConstantsHandle<ParticleDataTable>& pdt);

	 double LRUCorrection(double normalizedPosZ, double energy,int crystalId, ConditionsHandle<CalorimeterCalibrations>&);

//2016-03-18 Gianipez: maybe add also Poisson fluctuations?

	 double BirkesCorrection(int pdgId, double energy);



  };

  //--------------------------------------------------------------------

  void MakeCaloReadoutHitsDev2::beginJob()
  {
      art::ServiceHandle<art::TFileService> tfs;
      _hEtot      = tfs->make<TH1F>("hEtot",  "total energy deposited",        150,  0, 15);
      _hEner      = tfs->make<TH1F>("hEner",  "energy deposited / crystal",    200,  0, 20);
      _hinit      = tfs->make<TH1F>("hInit",  "Initial number of calo showers", 25,  0, 25);
      _hfinal     = tfs->make<TH1F>("hFinal", "Final number of calo showers",   25,  0, 25);
 }





  void  MakeCaloReadoutHitsDev2::produce(art::Event& event) 
  {
     art::ServiceHandle<GeometryService> geom;    
     if( !(geom->hasElement<Calorimeter>()) ) return;


     //update condition cache
     ConditionsHandle<AcceleratorParams> accPar("ignored");
     _mbtime = accPar->deBuncherPeriod;
     _toff.updateMap(event);


     // A container to hold the output hits.
     std::unique_ptr<CaloShowerCollection>  caloShowers(new CaloShowerCollection);


     //Get handles to caloStepMC collection
     art::Handle<CaloShowerStepMCCollection> caloShowerStepMCHandle,caloROShowerStepMCHandle;
     event.getByLabel(_caloShowerStepMCModuleLabel, _caloShowerMCName,   caloShowerStepMCHandle);


     makeCalorimeterHits(caloShowerStepMCHandle, *caloShowers);

     // Add the output hit collection to the event
     event.put(std::move(caloShowers));
  } 








//------------------------------------------------------------------------------------------------------------------------- 
  void MakeCaloReadoutHitsDev2::makeCalorimeterHits (art::Handle<CaloShowerStepMCCollection>const& caloShowerStepMCCollHandle, 
                                                     CaloShowerCollection &caloShowers) {
    
    
     GlobalConstantsHandle<ParticleDataTable>  pdt;
     ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

     Calorimeter const &cal = *(GeomHandle<Calorimeter>());
     double cryhalflength   = cal.caloGeomInfo().crystalHalfLength();


     CaloShowerStepMCCollection const &caloShowerSteps(*caloShowerStepMCCollHandle);
     CaloShowerStepMC           const *caloShowerStepMCBase = &caloShowerSteps.front();



     //-----------------------------------------------------------------------
     //store corrected energy deposits for each crystal hits
     
     std::map<int, std::vector<CORHit> >  hitMapCrystal;    
     
     
     for (auto istep = caloShowerSteps.begin(); istep !=caloShowerSteps.end(); ++istep)
     {        
     
	  CaloShowerStepMC const& step = *istep;
	  
	  size_t idx = (&step - caloShowerStepMCBase);
	  art::Ptr<CaloShowerStepMC> stepPtr = art::Ptr<CaloShowerStepMC>(caloShowerStepMCCollHandle,idx);

	  int crystalId    = step.volumeId();
	  int pdgId        = step.simParticle()->pdgId();
	  double edep_corr = step.energy();
	  double posZ      = step.position().z();

	  
	  //apply all the corrections on energy deposition -> scintillation light
	  if (_caloNonLinCorrection) 
	       edep_corr = nonLinearityCorrection(crystalId, pdgId, edep_corr, calorimeterCalibrations, pdt);

	  if (_caloLRUcorrection)
	       edep_corr = LRUCorrection(posZ/cryhalflength, edep_corr, crystalId, calorimeterCalibrations);

	  if (_caloBirkesCorrection)
	       edep_corr = BirkesCorrection(pdgId,edep_corr);
	  	  
	  
	  // time folding, see docdb-3425
	  double hitTimeUnfolded = _toff.totalTimeOffset(step.simParticle())+step.time();
	  double hitTime         = fmod(hitTimeUnfolded,_mbtime);

	  if (hitTime > _blindTime) hitMapCrystal[crystalId].emplace_back(CORHit(stepPtr,edep_corr,hitTime));
     }
     



     //Group hits within short time window
     for (auto &kv : hitMapCrystal )
     {
         
	 int crystalId = kv.first;
	 std::vector<CORHit> &hits = kv.second;
	 std::vector<CaloShower> caloShowerBuild;
	 
	 
	 std::sort(hits.begin(),hits.end(), [](CORHit const& a, CORHit const& b){return a.time_ < b.time_;});
	 
	 
	 double h_time  = hits[0].time_;
	 double h_edep  = hits[0].edep_corr_;
	 std::vector<art::Ptr<CaloShowerStepMC>> h_sims {hits[0].step_};
	   
	 for( size_t i=1; i<hits.size(); ++i )
	 {	  
	      //Save hit and create new hits
	      if ( (hits[i].time_ - h_time) > _timeGap )
	      {	    
	           if(_diagLevel > 1) std::cout<<"MakecaloReadoutHits2Dev:: insert hit crystalId="<<crystalId
		        	 	       <<" with time="<< h_time<<"and energy = "<<h_edep<<std::endl;		      		    
	  
	          _hEner->Fill(h_edep);

	           caloShowerBuild.push_back(CaloShower(crystalId,h_sims,h_time,h_edep));
	           h_sims.clear();
	           h_time    = hits[i].time_;
	           h_edep    = hits[i].edep_corr_;
		   h_sims.push_back(hits[i].step_);
	      } 	  
	      else
	      {
	           h_edep += hits[i].edep_corr_;
	           h_sims.push_back(hits[i].step_);
	      }
	 }
	 
	 //do not forget to save the last hit!
         caloShowerBuild.push_back(CaloShower(crystalId,h_sims,h_time,h_edep));
         if(_diagLevel > 1) std::cout<<"MakecaloReadoutHits2Dev:: insert hit crystalId="<<crystalId
	         	 	     <<" with time="<< h_time<<"and energy = "<<h_edep<<std::endl;		      		    	  
         _hinit->Fill(caloShowerBuild.size());

     
         //now a little bit of pre-filtering
	 
	 auto iterBefore  = caloShowerBuild.rend();
	 auto iterCurrent = caloShowerBuild.begin();
	 auto iterAfter   = caloShowerBuild.begin();

	 while(iterCurrent != caloShowerBuild.end())
	 {
	       ++iterAfter;     

	       double deltaTimePlus(1e6), deltaTimeMinus(1e6);
	       if (iterAfter   != caloShowerBuild.end() )  deltaTimePlus  = iterAfter->time() - iterCurrent->time();
	       if (iterBefore  != caloShowerBuild.rend())  deltaTimeMinus = iterCurrent->time() - iterBefore->time();

	       if (iterCurrent->energy() < _filterEnergy && std::min(deltaTimePlus,deltaTimeMinus) > _filterDeltaT)
	       {
	           if(_diagLevel > 1) std::cout<<" erase "<<crystalId<<" "<<iterCurrent->energy()<<" "<<iterCurrent->time()<<" "<<std::endl;
		   iterCurrent = caloShowerBuild.erase(iterCurrent);		      
	       } 
	       else
	       {
	           ++iterBefore;
		   ++iterCurrent;		
	       } 
	 }      

         _hfinal->Fill(caloShowerBuild.size());
         if (_diagLevel>1)  for (auto &cs : caloShowerBuild)  std::cout<<" final "<<crystalId<<" "<<cs.energy()<<" "<<cs.time()<<std::endl;
        
	
	 //finally copy the remaining caloshowers to the final vector 
         std::copy(caloShowerBuild.begin(),caloShowerBuild.end(),std::back_inserter(caloShowers));
     }

  }












  

 
  //-----------------------------------------------------------------------------
  double MakeCaloReadoutHitsDev2::nonLinearityCorrection(int cryId, int particleCode, double energy, 
                                                        ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations, 
						        GlobalConstantsHandle<ParticleDataTable>& pdt  )
  { 
    double edep(energy);
    //double MeV2keV = 1000.0;
      

    if (particleCode!= 11 && particleCode != 22 ) return edep;    
    /*    
    const HepPDT::ParticleData& data = pdt->particle(particleCode).ref();
    double mass = data.mass().value();
    double kinetic_energy = std::sqrt(h.momentum().mag2() + mass*mass) - mass;
    if (kinetic_energy > 1.0) return edep;
    
    //the formula used returns positive values if tmpEnergy>(approximately) 3keV, see Mu2e doc 1748-v1
    kinetic_energy = MeV2keV*kinetic_energy;
    double trackKine = kinetic_energy;
    kinetic_energy = std::log10(kinetic_energy);
    kinetic_energy *= _randGauss.fire(calorimeterCalibrations->LINpar2(cryId), calorimeterCalibrations->LINpar2Err(cryId) );
    kinetic_energy += _randGauss.fire(calorimeterCalibrations->LINpar1(cryId), calorimeterCalibrations->LINpar1Err(cryId) );
    kinetic_energy = std::log10(kinetic_energy);
    kinetic_energy *= _randGauss.fire(calorimeterCalibrations->LINpar0(cryId), calorimeterCalibrations->LINpar0Err(cryId) );
    kinetic_energy /= _randGauss.fire(calorimeterCalibrations->LINpar3(cryId), calorimeterCalibrations->LINpar3Err(cryId) );
    kinetic_energy /= std::log10(trackKine);
    
    if (kinetic_energy>0) energy *= kinetic_energy;	
    
    if (_diagLevel > 2) {
      std::cout<<"************************** BEFORE / AFTER NON-LINEARITY EFFECT-> edep_corr = "
	       << edep_save<<"  /  "<<energy<<std::endl
	       <<", energyKin = "
	       << trackKine 
	       << ", mass = "<< mass
	       <<std::endl;
    }
    */
    return edep;
    
  }

  //-----------------------------------------------------------------------------
  // apply a correction of type Energy = ((1-s)/L'*Z+s)*energy where Z position along the crystal, L' is the crystal half-length
  // and s is the intercept at Z=0 (i.e. non-uniformity facto, e.g. 5% -> s = 1.05)
  
  double MakeCaloReadoutHitsDev2::LRUCorrection(double normalizedPosZ, double energy, int crystalId, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
     double edep(energy);     
     
     double alpha  = calorimeterCalibrations->LRUpar0(crystalId);
     double factor =  (1.0+alpha) - alpha*normalizedPosZ;
     edep *= factor;
     if (_diagLevel > 2) std::cout<<"MakeCaloReadoutHiteDev2 before / after LRU -> edep_corr = "<< edep<<"  /  "<<energy<<"  at position Z="<<normalizedPosZ<<std::endl;	  
     return edep;     
  }



  //-----------------------------------------------------------------------------
  double MakeCaloReadoutHitsDev2::BirkesCorrection(int particleCode, double energy)
  {
     double edep(energy);
     if (particleCode==2212 || particleCode == 2112) edep/=4;
     return edep;    
  }




} 

using mu2e::MakeCaloReadoutHitsDev2;
DEFINE_ART_MODULE(MakeCaloReadoutHitsDev2);




