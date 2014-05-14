//
// An EDProducer Module that reads CaloHit objects and turns them into
// CaloCrystalHit objects, collection
//
// $Id: MakeCaloCrystalHitsNew_module.cc,v 1.7 2014/05/14 18:17:16 murat Exp $
// $Author: murat $
// $Date: 2014/05/14 18:17:16 $
//
// Original author KLG
// for realistic modeling: reduce timeGap from 30ns to 1 ns in Mu2eG4/test/calorimeter.txt
//
// Modified version
// Added pile-up management: two hits in the same crystal within the
// leading edge time are always merged, otherwise we do some geometrical 
// considerations (see code for details).
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
#include "ConditionsService/inc/GlobalConstantsHandle.hh"
#include "ConditionsService/inc/CalorimeterCalibrations.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "SeedService/inc/SeedService.hh"

#include "TH1F.h"
#include "TH2F.h"

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandGaussQ.h"


using namespace std;
using art::Event;


namespace mu2e {

  class MakeCaloCrystalHitsNew : public art::EDProducer {

  public:

    explicit MakeCaloCrystalHitsNew(fhicl::ParameterSet const& pset) :

      // Parameters

      _diagLevel                  (pset.get<int>   ("diagLevel"     ,0)),
      _drawLevel                  (pset.get<int>   ("drawLevel"     ,0)),
      _maxFullPrint               (pset.get<int>   ("maxFullPrint"  ,5)),
      _tDecay                     (pset.get<double>("tDecay"        , 40)),      // ns
      _tRise                      (pset.get<double>("tRise"         , 10)),       // ns
      _k                          (pset.get<double>("k"             , 1)), // k constant for the pile-up solving
      _minimumEnergy              (pset.get<double>("minimumEnergy" ,0.0001)), // MeV
      _maximumEnergy              (pset.get<double>("maximumEnergy" ,1000.0)), // MeV
      _minimumTimeGap             (pset.get<double>("minimumTimeGap",1.0)), // ns

      _caloChargeProductionEffects(pset.get<bool>  ("caloChargeProductionEffects",0)),
      _caloROnoiseEffect          (pset.get<bool>  ("caloROnoiseEffect"          ,0)),

      _g4ModuleLabel              (pset.get<string>     ("g4ModuleLabel")),
      _caloReadoutModuleLabel     (pset.get<std::string>("caloReadoutModuleLabel", "CaloReadoutHitsMaker")),

      _engine                     (createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
      _randPoisson                ( _engine ),
      _randGauss                  (_engine),
      _messageCategory            ("CaloHitMaker")

    {
      // Tell the framework what we make.
      produces<CaloCrystalHitCollection>();
    }
    virtual ~MakeCaloCrystalHitsNew() { }


    virtual void beginJob();
    void produce( art::Event& e);


  private:

    int _diagLevel;
    int _drawLevel;
    int _maxFullPrint;

    double _tDecay;
    double _tRise;
    double _k;
    double _tLeadingEdge;

    double _minimumEnergy;  // minimum energy in the RO to count it
    double _maximumEnergy;  // energy of a saturated RO
    double _minimumTimeGap; // to merge the hits
 
    bool   _caloChargeProductionEffects;
    bool   _caloROnoiseEffect;

    string _g4ModuleLabel;  // Name of the module that made the input hits.
    string _caloReadoutModuleLabel; // Name of the module that made the calo hits.
    CLHEP::HepRandomEngine& _engine;
    CLHEP::RandPoisson _randPoisson;
    CLHEP::RandGaussQ _randGauss;

    const std::string _messageCategory;

    TH1F* _hCaloHitTime;
    TH1F* _hDeltaT;
    TH1F* _hEnergy;
    TH1F* _hMergedDeltaT;
    TH1F* _hMergedAmplitudes; 
    TH2F* _hNotMergedAmplitudesDeltaT;
    TH2F* _hMergedAmplitudesDeltaT;
    
    void makeCrystalHits(CaloCrystalHitCollection& caloCrystalHits, art::Handle<CaloHitCollection>& caloHitsHandle) ;
    void fixEnergy(CaloCrystalHit& caloCrystalHit, int tnro, double electronEdep);
    void chargeProductionCorrection(CaloCrystalHit& caloCrystalHit,int roid, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);
    void readoutNoiseCorrection(CaloCrystalHit& caloCrystalHit, int roid, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations);
    double waveform(double A, double t0, double t);
  };

  void MakeCaloCrystalHitsNew::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    
    _hCaloHitTime      = tfs->make<TH1F>( "hCaloHitTime",  "Time of calorimeter hits; Time [ns]", 200, 0., 2000.);
    _hDeltaT           = tfs->make<TH1F>( "hDeltaT",  "Signals #Deltat; #Deltat [ns]", 200, 0., 2000.);
    _hEnergy           = tfs->make<TH1F>( "hEnergy",  "Energy of the signals; Energy [MeV]", 100, 0., 110);
    _hMergedDeltaT     = tfs->make<TH1F>( "hMergedDeltaT",  "Merged signals #DeltaT; #Deltat [ns]", 200, 0., 2000.);
    _hMergedAmplitudes = tfs->make<TH1F>( "hMergedAmplitudes",  "Ratio of amplitudes for merged signals; B/A", 200, 0., 1);
    _hNotMergedAmplitudesDeltaT = tfs->make<TH2F>("hNotMergedAmplitudesDeltaT", "B/A vs #Deltat for not merged signals", 200, 0, 1.2, 200, 0, 2000);      
    _hMergedAmplitudesDeltaT = tfs->make<TH2F>("hMergedAmplitudesDeltaT", "B/A vs #Deltat for merged signals", 200, 0, 1.2, 200, 0, 2000);
  }




  void MakeCaloCrystalHitsNew::produce(art::Event& event) {

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




//-----------------------------------------------------------------------------
  void MakeCaloCrystalHitsNew::makeCrystalHits(CaloCrystalHitCollection&        caloCrystalHits,
					       art::Handle<CaloHitCollection>&  caloHitsHandle ) {

    CaloHitCollection const& caloHits(*caloHitsHandle);
    if (caloHits.size()<=0) return;

    // Handle to the conditions servicea
    ConditionsHandle<CalorimeterCalibrations> calorimeterCalibrations("ignored");
    
     
    Calorimeter const & cal = *(GeomHandle<Calorimeter>());
    int nro                 = cal.nROPerCrystal();
    double electronEdep     = cal.getElectronEdep();


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
    for ( CaloHitCollection::const_iterator i=caloHits.begin(); i!=caloHits.end(); ++i ) {
      //cout << "Time of the hit " <<(*i).time() << endl;
      caloHitsSorted.push_back( &(*i));
    }
     
    sort ( caloHitsSorted.begin(), caloHitsSorted.end(), lessByCIdAndTimeByPointer<CaloHit>(&cal) );
     
     if ( _diagLevel > 2 ) {
      cout << __func__ << ": Total number of hit RO = " << caloHitsSorted.size() << endl;
      for( size_t i=0; i<caloHitsSorted.size(); ++i ) 
	cout << __func__ << ": " << caloHitsSorted[i]
	     << " Ro ID: " << caloHitsSorted[i]->id() 
	     << " CrystalId: " << cal.crystalByRO(caloHitsSorted[i]->id()) << endl;
    }

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
    } 
    else {
      caloCrystalHit.assignEnergyToTot(cal.crystalByRO(hit0.id()),hit0);
    }
     
    // Time of the first maximum in the two-signals function
    _tLeadingEdge = _tDecay*_tRise/(_tDecay-_tRise)*log(_tDecay/_tRise);
     
    CaloCrystalHitCollection signals;
    signals.reserve(caloHits.size());
    
    for( std::vector<CaloHit const *>::const_iterator i = caloHitsSorted.begin()+1; 
	 i != caloHitsSorted.end(); ++i) {
      CaloHit const& hit = **i;
      int cid = cal.crystalByRO(hit.id());
      
      if (_drawLevel) _hCaloHitTime->Fill(hit.time());
      
      if (_diagLevel) {
	cout << __func__ << ": Original RO hit: " << hit << endl;
	cout << __func__ << ": old, new cid:  " << caloCrystalHit.id() << ", " << cid << endl;
	cout << __func__ << ": old, new time: " << caloCrystalHit.time() << ", " << hit.time() << endl;
	cout << __func__ << ": time difference, gap: " << (hit.time() - caloCrystalHit.time()) << ", "<< _minimumTimeGap << endl;
      }
      
      // merge the hits if they are of the same crystal and closer than _tLeadingEdge
      if (caloCrystalHit.id() == cid && (( hit.time() - caloCrystalHit.time()) < _tLeadingEdge) ) {
	if ( hit.energyDep()>= _minimumEnergy && hit.energyDep() < _maximumEnergy ) {
	  size_t idx = ( &hit - base );
	  caloCrystalHit.add( hit, art::Ptr<CaloHit>(caloHitsHandle,idx));
	}
	else {
	  caloCrystalHit.addEnergyToTot(hit);
	}
	
	if (_diagLevel) cout << __func__ << ": Added to the hit:  " << caloCrystalHit << endl;
      } 
      else {
	if (caloCrystalHit.energyDep()>0.0) {
	  if (_diagLevel) cout << __func__ << ": Inserting old hit: " << caloCrystalHit << endl;
	  
	  int roId = hit.id();
	  if (_caloChargeProductionEffects) {
	    double esave = caloCrystalHit.energyDep();
	    chargeProductionCorrection(caloCrystalHit, roId, calorimeterCalibrations);
	    if (_diagLevel) { 
	      std::cout<<"before / after ChargeProductionEffects-correction, energy = "
		       << esave<<" / "<<caloCrystalHit.energyDep()<<std::endl;
	    }
	  }
	  
	  if (_caloROnoiseEffect) {
	    double esave = caloCrystalHit.energyDep();
	    readoutNoiseCorrection(caloCrystalHit,roId, calorimeterCalibrations);
	    if(_diagLevel) std::cout<<"after ROnoiseEffects-correction, energy = "
				    <<  esave<<" / "<<caloCrystalHit.energyDep()<<std::endl;	      
	  }
	  signals.push_back(caloCrystalHit);
	}
	   
	// this resets the caloCrystalHit and sets its id and puts one hit in
	if ( hit.energyDep()>= _minimumEnergy && hit.energyDep() < _maximumEnergy ) {
	  size_t idx = ( &hit - base );
	  caloCrystalHit.assign(cid, hit, art::Ptr<CaloHit>(caloHitsHandle,idx));
	} 
	else {
	  caloCrystalHit.assignEnergyToTot(cid,hit);
	}
	
	if (_diagLevel) cout << __func__ << ": Created new hit: " << caloCrystalHit << endl;
      }
    }
    
    if (caloCrystalHit.energyDep()>0.0) {
      if (_diagLevel) cout << __func__ << ": Inserting last old hit:" << caloCrystalHit << endl;
      signals.push_back(caloCrystalHit);
    }
    
    double   A, B;
    double   t_shift;
    double   t_second_peak;
    double   waveform_A, waveform_B, waveform_tot;
    
    CaloCrystalHit *firstSignal, *secondSignal;
    std::vector<CaloCrystalHit>::iterator j;
    
    if (_drawLevel) {
      _hNotMergedAmplitudesDeltaT->SetMarkerColor(kGreen);
      _hMergedAmplitudesDeltaT->SetMarkerColor(kRed);
      _hMergedAmplitudesDeltaT->GetXaxis()->SetTitle("B/A");
      _hMergedAmplitudesDeltaT->GetYaxis()->SetTitle("#Deltat [ns]");
      _hNotMergedAmplitudesDeltaT->GetXaxis()->SetTitle("B/A");
      _hNotMergedAmplitudesDeltaT->GetYaxis()->SetTitle("#Deltat [ns]");
    }
     
    for (std::vector<CaloCrystalHit>::iterator i = signals.begin(); i != signals.end(); ++i) {
      fixEnergy(*i,nro,electronEdep);
    }
     
    for (std::vector<CaloCrystalHit>::iterator i = signals.begin(); i != signals.end(); i=j) {
      if (i==signals.end()) break;
      
      firstSignal = i.operator->();
      j = i+1;
      
      if (j==signals.end()) {
	caloCrystalHits.push_back(*firstSignal); 
	break; 
      }
      
      secondSignal = j.operator->();
      
      if (_drawLevel) _hEnergy->Fill(firstSignal->energyDep());
      
      if (_diagLevel) {
	cout << "ID: " << firstSignal->id() << endl;
	cout << "Time: " << firstSignal->time() << endl;
      }
//-----------------------------------------------------------------------------
// while the second signal belongs to the same crystal of the first one
//-----------------------------------------------------------------------------
      while (firstSignal->id() == secondSignal->id()) {
	if (j == signals.end()) break;
	
	A       = firstSignal->energyDep()/(_tDecay-_tRise);
	B       = secondSignal->energyDep()/(_tDecay-_tRise);
	t_shift = secondSignal->time() - firstSignal->time();
	
	// Time of the second peak in the two-signals function
	t_second_peak =  - (log(_tDecay/_tRise) + t_shift/_tRise + 
			    log(B/A+exp(-t_shift/_tRise)) - t_shift/_tDecay - 
			    log(B/A+exp(-t_shift/_tDecay)))*_tDecay*_tRise/(_tRise-_tDecay);
	
	if (_diagLevel) {
	  cout << "SAME SIGNAL  : " << secondSignal->id() << endl;
	  cout << "T second peak: " << t_second_peak      << endl;
	}
	
	if (_drawLevel) {
	  _hEnergy->Fill(secondSignal->energyDep());
	  _hDeltaT->Fill(t_shift);
	}
//-----------------------------------------------------------------------------
// Values of the signals at the time of the second peak
//-----------------------------------------------------------------------------
	waveform_A   = waveform(A, 0, t_second_peak);
	waveform_B   = waveform(B, t_shift, t_second_peak);
	waveform_tot = waveform_A + waveform_B;
	
//-----------------------------------------------------------------------------
// check if pile-up is resolvable
//-----------------------------------------------------------------------------
	if ((waveform_tot - waveform_A) > _k*waveform_A) {
	  if (_drawLevel) _hNotMergedAmplitudesDeltaT->Fill(B/A, t_shift);
	  
	  caloCrystalHits.push_back(*firstSignal);
	  *firstSignal = *secondSignal;
	} 
	else {
	  if (_drawLevel) {
	    _hMergedAmplitudes->Fill(B/A);
	    _hMergedDeltaT->Fill(t_shift);
	    _hMergedAmplitudesDeltaT->Fill(B/A, t_shift);
	  }
	  
	  // add all the hits of the second signal to the first one (merging)

	  firstSignal->add(secondSignal);
	}
	
	if (j!=signals.end()) {
	  j++;
	  secondSignal = j.operator->();
	}
      }
      caloCrystalHits.push_back(*firstSignal);
    }    
  } 
  
  
  
  
  void MakeCaloCrystalHitsNew::chargeProductionCorrection(CaloCrystalHit& caloCrystalHit,int roId, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations){
     
     double energy     = caloCrystalHit.energyDep();
     double lightYield = _randGauss.fire(calorimeterCalibrations->ROpe(roId), calorimeterCalibrations->ROpeErr(roId));
     if (lightYield <= 0.0) lightYield = calorimeterCalibrations->ROpe(roId);
     
     double mean  = energy*lightYield/calorimeterCalibrations->ROfano(roId);
     double N_pe1 = _randPoisson.fire(mean);//gRandom->PoissonD(mean);//
     energy       = N_pe1*calorimeterCalibrations->ROfano(roId)/lightYield;//calorimeterCalibrations->APDpe(roId);  //Poissonian smearing: photostatistic

     caloCrystalHit.setEnergyDep(energy);
  }
  
  
  void MakeCaloCrystalHitsNew::readoutNoiseCorrection(CaloCrystalHit& caloCrystalHit,int roId, ConditionsHandle<CalorimeterCalibrations>& calorimeterCalibrations){
     
     double energy = caloCrystalHit.energyDep();
     double tmpS   = _randGauss.fire(0.0,  calorimeterCalibrations->ROnoise(roId) );
     energy += tmpS;
     if (energy > 0.0) caloCrystalHit.setEnergyDep(energy);
     else caloCrystalHit.setEnergyDep(0);
  }
  
   
  
  
  void MakeCaloCrystalHitsNew::fixEnergy(CaloCrystalHit& caloCrystalHit, int tnro, double electronEdep) {

    int nridu = caloCrystalHit.numberOfROIdsUsed();

    if ( _diagLevel > 0 ) {
      cout << __func__ << ": fixing energy: " << caloCrystalHit.energyDep()
	   << ", used roids: " << nridu << ", energyDepT: " << caloCrystalHit.energyDepTotal() 
	   << endl;
    }
//-----------------------------------------------------------------------------
// this was plain wrong - if N hits are merged into one, still need to divide by 
// the number of photosensors...
//-----------------------------------------------------------------------------
//      caloCrystalHit.setEnergyDep(caloCrystalHit.energyDep()/double(nridu));
//      caloCrystalHit.setEnergyDepTotal(caloCrystalHit.energyDepTotal()-float(nridu-1)*caloCrystalHit.energyDep());

    caloCrystalHit.setEnergyDep(caloCrystalHit.energyDep()/2.);
    caloCrystalHit.setEnergyDepTotal(caloCrystalHit.energyDepTotal()/2.);

    // fix only if all ro are saturated
    if (nridu == 0 && caloCrystalHit.energyDepTotal()/tnro >= electronEdep) {
      caloCrystalHit.setEnergyDep(caloCrystalHit.energyDepTotal());
    }
      
    if ( _diagLevel > 0 ) {
      cout << __func__ << ": fixed energy: " <<  caloCrystalHit.energyDep()
	   << ", used roids: " << nridu << ", energyDepT: " << caloCrystalHit.energyDepTotal() 
	   << endl;
    }

    return;
  }


  double MakeCaloCrystalHitsNew::waveform(double A, double t0, double t) {
    // Function which simulates the shape of the signal from the crystal
    return A*(exp(-(t-t0)/_tDecay)-exp(-(t-t0)/_tRise));
  }


}

using mu2e::MakeCaloCrystalHitsNew;
DEFINE_ART_MODULE(MakeCaloCrystalHitsNew);
