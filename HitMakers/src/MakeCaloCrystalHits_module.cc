//
// An EDProducer Module that reads CaloHit objects and turns them into
// CaloCrystalHit objects, collection
//
// $Id: MakeCaloCrystalHits_module.cc,v 1.17 2014/08/01 20:57:45 echenard Exp $
// $Author: echenard $
// $Date: 2014/08/01 20:57:45 $
//
// Original author KLG
//

// C++ includes.
#include <iostream>
#include <string>
#include <cmath>

// Framework includes.
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"

// Mu2e includes.
#include "CalorimeterGeom/inc/Calorimeter.hh"
#include "CalorimeterGeom/inc/sort_functors.hh"
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "SeedService/inc/SeedService.hh"


#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGaussQ.h"


using namespace std;
using art::Event;


namespace mu2e {

  class MakeCaloCrystalHits : public art::EDProducer {

  public:

    explicit MakeCaloCrystalHits(fhicl::ParameterSet const& pset) :

      // Parameters

      _diagLevel(pset.get<int>("diagLevel",0)),
      _maxFullPrint(pset.get<int>("maxFullPrint",5)),
      _minimumEnergy(pset.get<double>("minimumEnergy",0.0001)), // MeV
      _maximumEnergy(pset.get<double>("maximumEnergy",1000.0)), //MeV
      _minimumTimeGap(pset.get<double>("minimumTimeGap",1.0)),// ns
      _caloChargeProductionEffects(pset.get<bool>("caloChargeProductionEffects",0)),
      _caloCalibNoise(pset.get<double>("caloCalibNoise",0)),
      _caloROnoiseEffect(pset.get<bool>("caloROnoiseEffect",0)),
      _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
      _caloReadoutModuleLabel(pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),
      _engine( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
      _randPoisson( _engine ),
      _randGauss(_engine),
      _messageCategory("CaloHitMaker")

    {
      // Tell the framework what we make.
      produces<CaloCrystalHitCollection>();
    }
    virtual ~MakeCaloCrystalHits() { }


    virtual void beginJob();
    void produce( art::Event& e);


  private:

    int _diagLevel;
    int _maxFullPrint;

    double _minimumEnergy;  // minimum energy in the RO to count it
    double _maximumEnergy;  // energy of a saturated RO
    double _minimumTimeGap; // to merge the hits
    bool   _caloChargeProductionEffects;
    double _caloCalibNoise;
    bool   _caloROnoiseEffect;

    string _g4ModuleLabel;  // Name of the module that made the input hits.
    string _caloReadoutModuleLabel; // Name of the module that made the calo hits.
    CLHEP::HepRandomEngine& _engine;
    CLHEP::RandPoisson _randPoisson;
    CLHEP::RandGaussQ _randGauss;

    const std::string _messageCategory;

    void makeCrystalHits(CaloCrystalHitCollection& caloCrystalHits, art::Handle<CaloHitCollection>& caloHitsHandle) ;
    void fixEnergy(CaloCrystalHit& caloCrystalHit, int tnro, double electronEdep);
    void chargeProductionCorrection(CaloCrystalHit& caloCrystalHit,int roid, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);
    void readoutNoiseCorrection(CaloCrystalHit& caloCrystalHit, int roid, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);

  };


  

  void MakeCaloCrystalHits::beginJob(){}




  void MakeCaloCrystalHits::produce(art::Event& event) {


     if ( _diagLevel > 0 ) cout << __func__ << ": begin" << endl;

     static int ncalls(0);
     ++ncalls;

     // Check that calorimeter geometry description exists
     art::ServiceHandle<GeometryService> geom;    
     if( !(geom->hasElement<Calorimeter>()) ) return;


     //Get handles to calorimeter RO (aka APD) collection
     art::Handle<CaloHitCollection> caloHitsHandle;
     event.getByLabel(_caloReadoutModuleLabel, caloHitsHandle);
     if ( !caloHitsHandle.isValid()) return;
     
    //Create a new CaloCrystalHit collection and fill it
     unique_ptr<CaloCrystalHitCollection> caloCrystalHits(new CaloCrystalHitCollection);
     makeCrystalHits(*caloCrystalHits,caloHitsHandle);



     if ( _diagLevel > 0 ) {
        for (std::vector<CaloCrystalHit>::iterator i = (*caloCrystalHits).begin(); i != (*caloCrystalHits).end(); ++i) 
          std::cout<<"I have produced the CaloCrystalHit "<<(*i)<<std::endl;
        cout << __func__ << ": caloCrystalHits.size() "<< caloCrystalHits->size() << endl;
        cout << __func__ << ": ncalls " << ncalls << endl;
        cout << __func__ << ": end" << endl;
     }

     event.put(std::move(caloCrystalHits));

     return;
  }





  void MakeCaloCrystalHits::makeCrystalHits(CaloCrystalHitCollection& caloCrystalHits,
                                            art::Handle<CaloHitCollection>& caloHitsHandle) {

     CaloHitCollection const& caloHits(*caloHitsHandle);
     if (caloHits.size()<=0) return;

     // Handle to the conditions servicea
     ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
    
     
     Calorimeter const & cal = *(GeomHandle<Calorimeter>());
     int nro                 = cal.caloGeomInfo().nROPerCrystal();
     double electronEdep     = cal.caloGeomInfo().electronEdep();


     if ( _diagLevel > 2 ) {
       cout << __func__ << ": Total number of hit RO = " << caloHits.size() << endl;
       for( size_t i=0; i<caloHits.size(); ++i ) 
           cout << __func__ << ": " << caloHits[i]
        	<< " Ro ID: " << caloHits[i].id() 
        	<< " CrystalId: " << cal.crystalByRO(caloHits[i].id()) << endl;
     }



     // Sort hits by crystal id ( not readout id! ) and time
     // Need one level of indirection since objects in the event are const.
     std::vector<CaloHit const*> caloHitsSorted;
     caloHitsSorted.reserve(caloHits.size());
     for ( CaloHitCollection::const_iterator i=caloHits.begin(); i!=caloHits.end(); ++i )
       caloHitsSorted.push_back( &(*i));

     sort ( caloHitsSorted.begin(), caloHitsSorted.end(), lessByCIdAndTimeByPointer<CaloHit>(&cal) );




    // generate the CaloCrystalHits. First,
    // add all RO hits of the same crystal / time together (must clearly improve to deal with pileup)
    // if energy between energy_min and energy_max (i.e. adding signal from APDs)
    // energyDepTotal will include all the energy


    CaloHit const* base = &caloHits.front();
    CaloHit const& hit0 = **caloHitsSorted.begin();

    CaloCrystalHit caloCrystalHit;
    if ( hit0.energyDep()>= _minimumEnergy && hit0.energyDep() < _maximumEnergy ) {
       size_t idx = ( &hit0 - base );
       caloCrystalHit.assign(cal.crystalByRO(hit0.id()), hit0, art::Ptr<CaloHit>(caloHitsHandle,idx));
    } else {
       caloCrystalHit.assignEnergyToTot(cal.crystalByRO(hit0.id()),hit0);
    }



    for( std::vector<CaloHit const *>::const_iterator i = caloHitsSorted.begin()+1; i != caloHitsSorted.end(); ++i) 
    {

	CaloHit const& hit = **i;
	int cid = cal.crystalByRO(hit.id());

	if (_diagLevel) 
	{
           cout << __func__ << ": Original RO hit: " << hit << endl;
           cout << __func__ << ": old, new cid:  " << caloCrystalHit.id() << ", " << cid << endl;
           cout << __func__ << ": old, new time: " << caloCrystalHit.time() << ", " << hit.time() << endl;
           cout << __func__ << ": time difference, gap: " << (hit.time() - caloCrystalHit.time()) << ", "<< _minimumTimeGap << endl;
	}

//	if (caloCrystalHit.id() == cid && (( hit.time() - caloCrystalHit.time()) < _minimumTimeGap) ) {
	if (caloCrystalHit.id() == cid && (( hit.time() - caloCrystalHit.time()) < 0.1) ) 
	{

             if ( hit.energyDep()>= _minimumEnergy && hit.energyDep() < _maximumEnergy ) 
	     {
               size_t idx = ( &hit - base );
               caloCrystalHit.add( hit, art::Ptr<CaloHit>(caloHitsHandle,idx));
             } else {
               caloCrystalHit.addEnergyToTot(hit);
             }
             if (_diagLevel) cout << __func__ << ": Added to the hit:  " << caloCrystalHit << endl;


	} else {

             if (caloCrystalHit.energyDep()>0.0) 
	     {
        	  if (_diagLevel) cout << __func__ << ": Inserting old hit: " << caloCrystalHit << endl;

		  
		  int roId = hit.id();

		  if (_caloChargeProductionEffects) 
		  {
		      double esave = caloCrystalHit.energyDep();
		      chargeProductionCorrection(caloCrystalHit, roId, calorimeterCalibrations);
		      if(_diagLevel) std::cout<<"before / after ChargeProductionEffects-correction, energy = "<< esave<<" / "<<caloCrystalHit.energyDep()<<std::endl;
		  }

		  if (_caloROnoiseEffect)
		  {
        	    double esave = caloCrystalHit.energyDep();
		    readoutNoiseCorrection(caloCrystalHit,roId, calorimeterCalibrations);
		    if (_diagLevel) std::cout<<"after ROnoiseEffects-correction, energy = "<<  esave<<" / "<<caloCrystalHit.energyDep()<<std::endl;	      
		  }

		  if (_caloCalibNoise > 1e-6)
		  {
		    double esave = caloCrystalHit.energyDep();
		    double energy = caloCrystalHit.energyDep()*_randGauss.fire(1.0, _caloCalibNoise );
		    caloCrystalHit.setEnergyDep( (energy > 0 ? energy : 0.0));
		    if(_diagLevel) std::cout<<"after noiseCalib-correction, energy = "<<  esave<<" / "<<caloCrystalHit.energyDep()<<std::endl;	      
		  }
		  
		  caloCrystalHits.push_back(caloCrystalHit);
             }

             // this resets the caloCrystalHit and sets its id and puts one hit in
             if ( hit.energyDep()>= _minimumEnergy && hit.energyDep() < _maximumEnergy ) {
                 size_t idx = ( &hit - base );
                 caloCrystalHit.assign(cid, hit, art::Ptr<CaloHit>(caloHitsHandle,idx));
             } else {
                 caloCrystalHit.assignEnergyToTot( cid, hit);
             }
             if (_diagLevel) cout << __func__ << ": Created new hit:   " << caloCrystalHit << endl;

	}

    }



    if (caloCrystalHit.energyDep()>0.0) 
    {
      if (_diagLevel) cout << __func__ << ": Inserting last old hit: " << caloCrystalHit << endl;

      int roId = caloCrystalHit.id();

      if (_caloChargeProductionEffects) 
      {
	  double esave = caloCrystalHit.energyDep();
	  chargeProductionCorrection(caloCrystalHit, roId, calorimeterCalibrations);
	  if(_diagLevel) std::cout<<"before / after ChargeProductionEffects-correction, last id="<<caloCrystalHit.id()<<" energy = "<< esave<<" / "<<caloCrystalHit.energyDep()<<"   ratio   "<<esave/caloCrystalHit.energyDep()<<std::endl;
      }

      if (_caloROnoiseEffect)
      {
        double esave = caloCrystalHit.energyDep();
	readoutNoiseCorrection(caloCrystalHit,roId, calorimeterCalibrations);
	if(_diagLevel) std::cout<<"after ROnoiseEffects-correction, energy = "<<  esave<<" / "<<caloCrystalHit.energyDep()<<std::endl;	      
      }

      if (_caloCalibNoise > 1e-6)
      {
	double esave = caloCrystalHit.energyDep();
	double energy = caloCrystalHit.energyDep()*_randGauss.fire(1.0, _caloCalibNoise );
	caloCrystalHit.setEnergyDep( (energy > 0 ? energy : 0.0));
	if(_diagLevel) std::cout<<"after noiseCalib-correction, energy = "<<  esave<<" / "<<caloCrystalHit.energyDep()<<std::endl;	      
      }


      caloCrystalHits.push_back(caloCrystalHit);
    }


    //sort hits by time, important for clustering
    std::sort(caloCrystalHits.begin(),caloCrystalHits.end(),[](CaloCrystalHit const& a, CaloCrystalHit const& b) {return a.time() < b.time();});	    


    //Each CaloCrystalHit has the sum of the signal APD, must rescale for the total ebergy
    for (std::vector<CaloCrystalHit>::iterator i = caloCrystalHits.begin(); i != caloCrystalHits.end(); ++i) 
      fixEnergy(*i,nro,electronEdep);


  } 
  
  
  
  
  void MakeCaloCrystalHits::chargeProductionCorrection(CaloCrystalHit& caloCrystalHit,int roId, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations){
     
     double energy     = caloCrystalHit.energyDep();
     double lightYield = _randGauss.fire(calorimeterCalibrations->ROpe(roId), calorimeterCalibrations->ROpeErr(roId));
     if (lightYield <= 0.0) lightYield = calorimeterCalibrations->ROpe(roId);
     
     double mean  = energy*lightYield/calorimeterCalibrations->ROfano(roId);
     double N_pe1 = _randPoisson.fire(mean);//gRandom->PoissonD(mean);//
     energy       = N_pe1*calorimeterCalibrations->ROfano(roId)/lightYield;//calorimeterCalibrations->APDpe(roId);  //Poissonian smearing: photostatistic

     caloCrystalHit.setEnergyDep(energy);
  }
  
  
  void MakeCaloCrystalHits::readoutNoiseCorrection(CaloCrystalHit& caloCrystalHit,int roId, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations){
     
     double energy = caloCrystalHit.energyDep();
     double tmpS   = _randGauss.fire(0.0,  calorimeterCalibrations->ROnoise(roId) );
     energy += tmpS;
     if (energy > 0.0) caloCrystalHit.setEnergyDep(energy);
     else caloCrystalHit.setEnergyDep(0);
  }
  
   
  
  
  void MakeCaloCrystalHits::fixEnergy(CaloCrystalHit& caloCrystalHit, int tnro, double electronEdep) {

      int nridu = caloCrystalHit.numberOfROIdsUsed();

      if ( _diagLevel > 0 ) 
	cout << __func__ << ": fixing energy: " << caloCrystalHit.energyDep()
             << ", used roids: " << nridu << ", energyDepT: " << caloCrystalHit.energyDepTotal() << endl;

      caloCrystalHit.setEnergyDep(caloCrystalHit.energyDep()/double(nridu));
      caloCrystalHit.setEnergyDepTotal(caloCrystalHit.energyDepTotal()-float(nridu-1)*caloCrystalHit.energyDep());

      // fix only if all ro are saturated
      if (nridu == 0 && caloCrystalHit.energyDepTotal()/tnro >= electronEdep) caloCrystalHit.setEnergyDep(caloCrystalHit.energyDepTotal());
      
      if ( _diagLevel > 0 ) 
	cout << __func__ << ": fixed  energy: " <<  caloCrystalHit.energyDep()
             << ", used roids: " << nridu << ", energyDepT: " << caloCrystalHit.energyDepTotal() << endl;
      

      return;

  }


}

using mu2e::MakeCaloCrystalHits;
DEFINE_ART_MODULE(MakeCaloCrystalHits);
