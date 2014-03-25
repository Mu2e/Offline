//
// This module transforms StrawDigi objects into StrawHit objects
// It also builds the truth match map (if MC truth info for the StrawDigis exists)
//
// $Id: StrawHitsFromStrawDigis_module.cc,v 1.12 2014/03/25 22:14:39 brownd Exp $
// $Author: brownd $ 
// $Date: 2014/03/25 22:14:39 $
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
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "TrackerConditions/inc/StrawElectronics.hh"
#include "TrackerConditions/inc/StrawPhysics.hh"
//CLHEP
// data
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"

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

    // # of ADC digitizations to sum to define baseline
    unsigned _nbase;
    double _mbtime; // period of 1 microbunch
    double _mbbuffer; // buffer on that for ghost hits (wrapping)
    double _maxdt; // maximum time difference between end times
    bool _singledigi; // turn single-end digitizations into hits
// Diagnostics level.
    int _printLevel, _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the StrawDigi collection
    std::string _strawDigis;
    ConditionsHandle<StrawElectronics> _strawele; // models of straw response to stimuli
    ConditionsHandle<StrawPhysics> _strawphys; // models of straw response to stimuli
  };

  StrawHitsFromStrawDigis::StrawHitsFromStrawDigis(fhicl::ParameterSet const& pset) :
    _nbase(pset.get<unsigned>("NumADCBaseline",1)),
    _mbbuffer(pset.get<double>("TimeBuffer",100.0)), // nsec
    _maxdt(pset.get<double>("MaxTimeDifference",8.0)), // nsec
    _singledigi(pset.get<bool>("UseSingleDigis",false)), // use or not single-end digitizations
    _printLevel(pset.get<int>("printLevel",0)),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _strawDigis(pset.get<string>("StrawDigis","makeSD"))
    {
      produces<StrawHitCollection>();
      produces<PtrStepPointMCVectorCollection>("StrawHitMCPtr");
      produces<StrawDigiMCCollection>("StrawHitMC");
      if(_printLevel > 0) cout << "In StrawHitsFromStrawDigis constructor " << endl;
    }

  void StrawHitsFromStrawDigis::beginJob(){
  }

  void StrawHitsFromStrawDigis::beginRun( art::Run& run ){
  }

  void StrawHitsFromStrawDigis::produce(art::Event& event) {
    if(_printLevel > 0) cout << "In StrawHitsFromStrawDigis produce " << endl;
// update conditions
    
    _strawele = ConditionsHandle<StrawElectronics>("ignored");
    _strawphys = ConditionsHandle<StrawPhysics>("ignored");
    unique_ptr<StrawHitCollection>             strawHits(new StrawHitCollection);
    unique_ptr<PtrStepPointMCVectorCollection> mcptrHits(new PtrStepPointMCVectorCollection);
    unique_ptr<StrawDigiMCCollection> mchits(new StrawDigiMCCollection);
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime = accPar->deBuncherPeriod;

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
    const StrawDigiMCCollection * mcdigis(0);
    art::Handle<StrawDigiMCCollection> mcdigiH;
    if(event.getByLabel(_strawDigis,"StrawDigiMC",mcdigiH))
      mcdigis = mcdigiH.product();
  // loop over digis.  Note the MC truth is in sequence
    size_t ndigi = strawdigis->size();
    if( (mcptrdigis != 0 && mcptrdigis->size() != ndigi) ||
	(mcdigis != 0 && mcdigis->size() != ndigi) )
      throw cet::exception("RECO")<<"mu2e::StrawHitsFromStrawDigis: MCPtrDigi collection size doesn't match StrawDigi collection size" << endl;
    for(size_t isd=0;isd<ndigi;++isd){
      StrawDigi const& digi = (*strawdigis)[isd];
// convert the digi to a hit
      array<double,2> times;
      _strawele->tdcTimes(digi.TDC(),times);
// hit wants primary time and dt.  Check if both ends digitized, or if
// this is a single-end digitization
      double time(times[0]);
      double dt = times[1]-times[0];
      bool makehit(true);
      if(time < _mbtime+_mbbuffer && fabs(dt)<_maxdt ){
	time = times[0];
      } else if(_singledigi){
// single-ended hit.  Take the valid time, and set delta_t to 0.  This needs
// to be flaged in StrawHit, FIXME!!!
	if(times[0] < _mbtime+_mbbuffer)
	  time = times[0];
	else if(times[1] < _mbtime+_mbbuffer)
	  time = times[1];
	else
	  makehit = false;
      } else
	makehit = false;
      if(makehit){
// to get the charge we should fit the whole waveform: for now, just integrage  minus the baseline
// FIXME!!
	static const double pcfactor(1.0e-3);
	StrawDigi::ADCWaveform const& adc = digi.adcWaveform();
	// note: pedestal is being subtracting inside strawele, in the real experiment we will need
	// per-channel version of this FIXME!!!
	double baseline(0.0);
	for(unsigned ibase=0;ibase<_nbase;++ibase){
	  baseline += _strawele->adcCurrent(adc[ibase]);
	}
	baseline /= double(_nbase);	
	double charge(0.0);
	for(size_t iadc=_nbase;iadc<adc.size();++iadc){
	  charge += (_strawele->adcCurrent(adc[iadc])-baseline)*_strawele->adcPeriod()*pcfactor;
	}
	// double the energy, since we only digitize 1 end of the straw.
	// use time division to correct for attenuation FIXME!!
	// the gain should come from a straw-dependent database FIXME!!
	double energy = 2.0*_strawphys->ionizationEnergy(charge/_strawphys->strawGain(2.0,0.0));
	// crate the straw hit and append it to the list
	StrawHit newhit(digi.strawIndex(),time,dt,energy);
	strawHits->push_back(newhit);
// copy MC truth from digi to hit.  These are exactly the same as for the digi
	if(mcptrdigis != 0){
	  mcptrHits->push_back((*mcptrdigis)[isd]);
	}
	if(mcdigis != 0){
	  mchits->push_back((*mcdigis)[isd]);
	}
      }
    }
// put objects into event
    event.put(std::move(strawHits));
    if(mcptrdigis != 0)event.put(std::move(mcptrHits),"StrawHitMCPtr");
    if(mchits != 0)event.put(move(mchits),"StrawHitMC");
  }
}
using mu2e::StrawHitsFromStrawDigis;
DEFINE_ART_MODULE(StrawHitsFromStrawDigis);

