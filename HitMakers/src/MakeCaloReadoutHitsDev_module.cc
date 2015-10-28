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
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/ParticleDataTable.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/CaloCrystalOnlyHitCollection.hh"
#include "MCDataProducts/inc/PDGCode.hh"
#include "MCDataProducts/inc/CaloShowerStepMCCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepMC.hh"
#include "Mu2eUtilities/inc/SimParticleTimeOffset.hh"
#include "SeedService/inc/SeedService.hh"

// Other includes.
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "TH1D.h"




namespace mu2e {



    // Anonymous namespace to hold some helper classes. COR is a compression of crystal or Readout
    namespace {

	 class CORHit {

	     public:
               
	       enum step_type {crystal, readout};

	       CORHit(art::Ptr<CaloShowerStepMC> const& step, step_type type, double edep_corr, double time) : 
	         _step(step), _type(type), _edep_corr(edep_corr), _time(time)
	       {}

	       art::Ptr<CaloShowerStepMC> _step;
	       step_type                  _type;
	       double                     _edep_corr;           
	       double                     _time;                
	 }; 
    } 








  //--------------------------------------------------------------------
  class MakeCaloReadoutHitsDev : public art::EDProducer {
 
    public:

	explicit MakeCaloReadoutHitsDev(fhicl::ParameterSet const& pset) :

	  // Parameters
	  _caloShowerStepMCModuleLabel(pset.get<std::string>("caloShowerStepMCModuleLabel")),
	  _caloShowerMCName(           pset.get<std::string>("caloShowerMCName")),
	  _caloROShowerMCName(         pset.get<std::string>("caloROShowerMCName")),
	  _caloLRUcorrection(          pset.get<bool>("caloLRUcorrection",0)),
	  _caloNonLinCorrection(       pset.get<bool>("caloNonLinCorrection",0)),
	  _caloBirkesCorrection(       pset.get<bool>("caloBirkesCorrection",0)),
	  _toff(                       pset.get<fhicl::ParameterSet>("TimeOffsets", fhicl::ParameterSet())),
	  _blindTime(                  pset.get<double>("blindTime",300.)),
	  _timeGap(                    pset.get<double>("timeGap",60.)),
	  _messageCategory(            "CaloReadoutHitsMakerDev"),
	  _diagLevel(                  pset.get<int>("diagLevel",0)),
	  _maxFullPrint(               pset.get<int>("maxFullPrint",5)),
	  _randGauss(                  createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
	  _hEtot(0),_hEner(0),_htime(0),_hZ(0)
	{
	  produces<CaloHitCollection>();
	  produces<CaloHitMCTruthCollection>();
	}

	virtual ~MakeCaloReadoutHitsDev() { }
	virtual void beginJob();
	void produce( art::Event& e);



     private:

	 std::string            _caloShowerStepMCModuleLabel;   
	 std::string            _caloShowerMCName;   
	 std::string            _caloROShowerMCName;   

	 bool                   _caloLRUcorrection;
	 bool                   _caloNonLinCorrection;
	 bool                   _caloBirkesCorrection;

	 SimParticleTimeOffset  _toff; 
	 double                 _blindTime;
	 double                 _mbtime;
	 double                 _timeGap;

	 const std::string      _messageCategory;
	 int                    _diagLevel;  
	 int                    _maxFullPrint;
	 CLHEP::RandGaussQ      _randGauss;

	 TH1F*                  _hEtot;
	 TH1F*                  _hEner;
	 TH1F*                  _htime;
	 TH1F*                  _hZ;




	 void makeCalorimeterHits (art::Handle<CaloShowerStepMCCollection>const&, 
                        	   art::Handle<CaloShowerStepMCCollection>const&, 
                        	   CaloHitCollection &,
				   CaloHitMCTruthCollection&);

	 double nonLinearityCorrection( int cryId, int particleCode, double energy, 
				      ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations,
				      GlobalConstantsHandle<ParticleDataTable>& pdt);

	 double LRUCorrection(double normalizedPosZ, double energy,int crid, ConditionsHandle<CalorimeterCalibrations>&);

	 double BirkesCorrection(int pdgId, double energy);



  };

  //--------------------------------------------------------------------

  void MakeCaloReadoutHitsDev::beginJob()
  {
      art::ServiceHandle<art::TFileService> tfs;
      _hEtot      = tfs->make<TH1F>("hEtot", "total energy deposited",      150,  0, 150);
      _hEner      = tfs->make<TH1F>("hEner", "energy deposiued / crystal", 2000,  0, 2000);
      _htime      = tfs->make<TH1F>("htime", "time of energy deposit",     2000,  0, 2000);
      _hZ         = tfs->make<TH1F>("hZ",    "Z",                           250,  0, 250);
  }





  void  MakeCaloReadoutHitsDev::produce(art::Event& event) 
  {
     art::ServiceHandle<GeometryService> geom;    
     if( !(geom->hasElement<Calorimeter>()) ) return;


     //update condition cache
     ConditionsHandle<AcceleratorParams> accPar("ignored");
     _mbtime = accPar->deBuncherPeriod;
     _toff.updateMap(event);


     // A container to hold the output hits.
     std::unique_ptr<CaloHitCollection>        caloHits(new CaloHitCollection);
     std::unique_ptr<CaloHitMCTruthCollection> caloMCHits(new CaloHitMCTruthCollection);


     //Get handles to caloStepMC collection
     art::Handle<CaloShowerStepMCCollection> caloShowerStepMCHandle,caloROShowerStepMCHandle;
     event.getByLabel(_caloShowerStepMCModuleLabel, _caloShowerMCName,   caloShowerStepMCHandle);
     event.getByLabel(_caloShowerStepMCModuleLabel, _caloROShowerMCName, caloROShowerStepMCHandle);


     makeCalorimeterHits(caloShowerStepMCHandle, caloROShowerStepMCHandle, *caloHits, *caloMCHits);

     // Add the output hit collection to the event
     event.put(std::move(caloHits));
     event.put(std::move(caloMCHits));
  } 








//------------------------------------------------------------------------------------------------------------------------- 
  void MakeCaloReadoutHitsDev::makeCalorimeterHits (art::Handle<CaloShowerStepMCCollection>const& caloShowerStepMCCollHandle, 
                                                    art::Handle<CaloShowerStepMCCollection>const& caloROShowerStepMCCollHandle,
                                                    CaloHitCollection &caloHits,
						    CaloHitMCTruthCollection& caloHitsMCTruth) {
    
    
     GlobalConstantsHandle<ParticleDataTable>  pdt;
     ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");

     Calorimeter const & cal = *(GeomHandle<Calorimeter>());
     double cryhalflength    = cal.caloGeomInfo().crystalHalfLength();



     CaloShowerStepMCCollection const& caloShowerSteps(*caloShowerStepMCCollHandle);
     CaloShowerStepMCCollection const& caloROShowerSteps(*caloROShowerStepMCCollHandle);

     CaloShowerStepMC           const* caloShowerStepMCBase   = &caloShowerSteps.front();
     CaloShowerStepMC           const* caloROShowerStepMCBase = &caloROShowerSteps.front();

     std::map<int, std::vector<CORHit> >  hitMap;    
     
     double totEner(0);
     for (auto istep = caloShowerSteps.begin(); istep !=caloShowerSteps.end(); ++istep)
     {        

	  CaloShowerStepMC const& step = *istep;
	  
	  size_t idx = (&step - caloShowerStepMCBase);
	  art::Ptr<CaloShowerStepMC> stepPtr = art::Ptr<CaloShowerStepMC>(caloShowerStepMCCollHandle,idx);

	  int crid         = step.volumeId();
	  int pdgId        = step.simParticle()->pdgId();
	  double edep_corr = step.energy();
	  double posZ      = step.position().z();

	  
	  //apply all the corrections on energy deposition -> scintillation light
	  if (_caloNonLinCorrection) 
	       edep_corr = nonLinearityCorrection(crid, pdgId, edep_corr, calorimeterCalibrations, pdt);

	  if (_caloLRUcorrection)
	       edep_corr = LRUCorrection(posZ/cryhalflength, edep_corr, crid, calorimeterCalibrations);

	  if (_caloBirkesCorrection)
	       edep_corr = BirkesCorrection(pdgId,edep_corr);
	  
	  // time folding, see docdb-3425
	  double hitTimeUnfolded = _toff.totalTimeOffset(step.simParticle())+step.time();
	  double hitTime         = fmod(hitTimeUnfolded,_mbtime);

	  if (hitTime > _blindTime)
	  {
              int ROidBase = cal.ROBaseByCrystal(crid);
              int ROidEnd  = ROidBase+cal.caloGeomInfo().nROPerCrystal();
	      for (int roid=ROidBase;roid<ROidEnd;++roid) hitMap[roid].push_back(CORHit(stepPtr, CORHit::crystal,edep_corr,hitTime));
	  }
	  
	  _hEner->Fill(crid,step.energy());
	  _htime->Fill(step.time(),step.nCompress());
	  _hZ->Fill(posZ,step.nCompress());
	  totEner +=  step.energy();   
     }
     _hEtot->Fill(totEner);
     

     for (auto istep = caloROShowerSteps.begin(); istep !=caloROShowerSteps.end(); ++istep)
     {        
	  CaloShowerStepMC const& step = *istep;
	  size_t idx = (&step - caloROShowerStepMCBase);
	  art::Ptr<CaloShowerStepMC> stepPtr = art::Ptr<CaloShowerStepMC>(caloROShowerStepMCCollHandle,idx);

	  int roid  = step.volumeId();

	  // time folding 
	  double hitTimeUnfolded = _toff.totalTimeOffset(step.simParticle())+step.time();
	  double hitTime         = fmod(hitTimeUnfolded,_mbtime);

	  if (hitTime > _blindTime) hitMap[roid].push_back(CORHit(stepPtr,CORHit::readout,step.energy(),hitTime));
     }




     for (auto& roIter : hitMap)
     {     
          int roid = roIter.first;
	  std::vector<CORHit> &hits = roIter.second;

          std::sort(hits.begin(), hits.end(),[](CORHit const& a, CORHit const& b) {return a._time < b._time;});


	  std::vector<art::Ptr<CaloShowerStepMC> > steps;
	  steps.push_back(hits[0]._step);
	  double            h_time = hits[0]._time;
	  double            h_edep = hits[0]._edep_corr;
	  CORHit::step_type h_type = hits[0]._type;


	  for( size_t i=1; i<hits.size(); ++i )
	  {
	       if( (hits[i]._time - h_time) > _timeGap )
	       {
		   //flush previous hit ...
		   if(_diagLevel > 1) std::cout<<"MakeCaloReadoutHitsDev:: insert hit roid="<<roid<<" with time="<< h_time
					       <<"and energy = "<<h_edep<<std::endl;		      		    

		   caloHits.push_back(       CaloHit(       roid,h_time,h_edep));
		   caloHitsMCTruth.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_type));

		   // ...and create new hit
		   steps.clear();
		   steps.push_back(hits[0]._step);

		   h_time    = hits[i]._time;
		   h_edep    = hits[i]._edep_corr;
		   h_type    = hits[i]._type;	    
	       } 
	       else
	       {
	          //just add the hit
   	          h_edep  += hits[i]._edep_corr;
                  if( hits[i]._type == CORHit::readout ) h_type = CORHit::readout; 
 	          steps.push_back(hits[i]._step);	      
	       }   
	  }

          //final flush
 	  if(_diagLevel > 1) std::cout<<"MakeCaloReadoutHitsDev:: insert hit roid="<<roid<<" with time="<< h_time
				      <<"and energy = "<<h_edep<<std::endl;		      		    

	  caloHits.push_back(       CaloHit(       roid,h_time,h_edep));
	  caloHitsMCTruth.push_back(CaloHitMCTruth(roid,h_time,h_edep,h_type));    
     }


  }












  

 
  //-----------------------------------------------------------------------------
  double MakeCaloReadoutHitsDev::nonLinearityCorrection(int cryId, int particleCode, double energy, 
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
  double MakeCaloReadoutHitsDev::LRUCorrection(double normalizedPosZ, double energy, int crid, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations)
  {
     double edep(energy);     
     
     double alpha  = calorimeterCalibrations->LRUpar0(crid);
     double factor =  (1.0+alpha) - alpha*normalizedPosZ;
     edep *= factor;
     if (_diagLevel > 2) std::cout<<"***************BEFORE / AFTER LRU EFFECT-> edep_corr = "<< edep<<"  /  "<<energy<<"  at position Z="<<normalizedPosZ<<std::endl;	  
     return edep;     
  }



  //-----------------------------------------------------------------------------
  double MakeCaloReadoutHitsDev::BirkesCorrection(int particleCode, double energy)
  {
     double edep(energy);
     if (particleCode==2212 || particleCode == 2112) edep/=4;
     return edep;    
  }




} 

using mu2e::MakeCaloReadoutHitsDev;
DEFINE_ART_MODULE(MakeCaloReadoutHitsDev);



