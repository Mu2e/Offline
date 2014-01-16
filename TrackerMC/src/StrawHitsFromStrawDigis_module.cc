//
// This module transforms StrawDigi objects into StrawHit objects
// It also builds the truth match map (if MC truth info for the StrawDigis exists)
//
// $Id: StrawHitsFromStrawDigis_module.cc,v 1.5 2014/01/16 21:03:22 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/01/16 21:03:22 $
//
// Original author David Brown, LBNL
//
// framework
#include "art/Framework/Principal/Event.h"
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Handle.h"
#include "GeometryService/inc/GeomHandle.hh"
#include "art/Framework/Core/EDProducer.h"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Services/Optional/TFileService.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "TrackerMC/inc/StrawElectronics.hh"
//CLHEP
#include "CLHEP/Random/RandGaussQ.h"
// data
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"

using namespace std;
namespace mu2e {

  class StrawHitsFromStrawDigis : public art::EDProducer {

  public:
    explicit StrawHitsFromStrawDigis(fhicl::ParameterSet const& pset);
    // Accept compiler written d'tor.

    virtual void beginJob();
    virtual void beginRun( art::Run& run );
    virtual void produce( art::Event& e);

  private:

    // Diagnostics level.
    int _printLevel, _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the StrawDigi collection
    std::string _strawDigis;
    double _EIonize; // Geant energy of each ionization (MeV)
    double _QIonize; // charge of a single ionization (=e, pC)
    double _gasgain; // avalanche gain
    double _G4dEdQ; // G4 equivalence of energy (MeV) and ionization charge (pCoulomb).  Should
    // come from a materials database, it depends FIXME!!
    StrawElectronics _strawele; // models of straw response to stimuli
  };

  StrawHitsFromStrawDigis::StrawHitsFromStrawDigis(fhicl::ParameterSet const& pset) :
    _printLevel(pset.get<int>("printLevel",1)),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _strawDigis(pset.get<string>("StrawDigis","makeSD")),
    _EIonize(pset.get<double>("EnergyPerIonization",100.0e-6)), // 27 ev/ionization for 100% Ar! , should use Ar/CO2 FIXME!!
    _QIonize(pset.get<double>("ChargePerIonization",1.6e-7)), // e, pC
    _gasgain(pset.get<double>("GasGain",3.0e4)),
    _strawele(pset.get<fhicl::ParameterSet>("StrawElectronics",fhicl::ParameterSet()))
    {
    // this should be computed per straw, including calibration effects FIXME!!!
      _G4dEdQ = _EIonize/(_QIonize*_gasgain);
      produces<StrawHitCollection>();
      produces<PtrStepPointMCVectorCollection>("StrawHitMCPtr");
      if(_printLevel > 0) cout << "In StrawHitsFromStrawDigis constructor " << endl;
    }

  void StrawHitsFromStrawDigis::beginJob(){
  }

  void StrawHitsFromStrawDigis::beginRun( art::Run& run ){
  }

  void StrawHitsFromStrawDigis::produce(art::Event& event) {
    if(_printLevel > 0) cout << "In StrawHitsFromStrawDigis produce " << endl;
    unique_ptr<StrawHitCollection>             strawHits(new StrawHitCollection);
    unique_ptr<PtrStepPointMCVectorCollection> mcptrHits(new PtrStepPointMCVectorCollection);

    // find the digis
    art::Handle<mu2e::StrawDigiCollection> strawdigisH; 
    const StrawDigiCollection* strawdigis(0);
    if(event.getByLabel(_strawDigis,strawdigisH))
      strawdigis = strawdigisH.product();
    if(strawdigis == 0)
      throw cet::exception("RECO")<<"mu2e::StrawHitsFromStrawDigis: No StrawDigi collection found for label " <<  _strawDigis << endl;

  // find the associated MC truth collection.  Note this doesn't have to exist!
    const PtrStepPointMCVectorCollection * mcptrdigis(0);
    art::Handle<PtrStepPointMCVectorCollection> mcptrdigiH;
    if(event.getByLabel(_strawDigis,"StrawDigiMCPtr",mcptrdigiH))
      mcptrdigis = mcptrdigiH.product();
  // loop over digis.  Note the MC truth is in sequence
    size_t ndigi = strawdigis->size();
    if(mcptrdigis != 0 && mcptrdigis->size() != ndigi)
      throw cet::exception("RECO")<<"mu2e::StrawHitsFromStrawDigis: MCPtrDigi collection size doesn't match StrawDigi collection size" << endl;
    for(size_t isd=0;isd<ndigi;++isd){
      StrawDigi const& digi = (*strawdigis)[isd];
// convert the digi to a hit
      array<double,2> times;
      _strawele.tdcTimes(digi.TDC(),times);
// hit wants primary time and dt.  Note we may have conflicting definitions of the sign of dt, FIXME!!
      double time = times[0];
      double dt = times[1]-times[0];
// to get the charge we should fit the whole waveform: for now, just integrage  minus the baseline
// FIXME!!
      static const double pcfactor(1.0e-3);
      StrawDigi::ADCWaveform const& adc = digi.adcWaveform();
      // note: pedestal is being subtracting inside strawele, in the real experiment we will need
      // per-channel version of this FIXME!!!
      double baseline = (_strawele.adcCurrent(adc[0]) + _strawele.adcCurrent(adc[1]))/2.0;
      double charge(0.0);
      for(size_t iadc=2;iadc<adc.size();++iadc){
	charge += (_strawele.adcCurrent(adc[iadc])-baseline)*_strawele.adcPeriod()*pcfactor;
      }
      // double the energy, since we only digitize 1 end of the straw
      double energy = 2.0*charge*_G4dEdQ;
  // crate the straw hit and append it to the list
      StrawHit newhit(digi.strawIndex(),time,dt,energy);
      strawHits->push_back(newhit);
// copy MC truth from digi to hit.  These are exactly the same as for the digi
      if(mcptrdigis != 0){
	mcptrHits->push_back((*mcptrdigis)[isd]);
      }
    }
// put objects into event
    event.put(std::move(strawHits));
    if(mcptrdigis != 0)event.put(std::move(mcptrHits),"StrawHitMCPtr");

  }
}
using mu2e::StrawHitsFromStrawDigis;
DEFINE_ART_MODULE(StrawHitsFromStrawDigis);

