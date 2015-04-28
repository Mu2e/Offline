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
// helpers
#include "TrkChargeReco/inc/PeakFit.hh"
#include "TrkChargeReco/inc/PeakFitRoot.hh"
#include "TrkChargeReco/inc/PeakFitFunction.hh"
//CLHEP
// root 
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TMarker.h"
#include "TTree.h"
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
    bool _truncateADC; // model ADC truncation
    bool _floatPedestal; // float pedestal in fit
    bool _floatWidth; // _float width in fit
    bool _earlyPeak, _latePeak, _findPeaks; // additional peak finding
// Diagnostics level.
    int _printLevel, _diagLevel, _debugLevel;
// Diagnostics
    TTree* _shdiag;
    Int_t _shdevice, _shsector, _shlayer, _shstraw, _shqstat;
    Float_t _shq, _shqt, _shqchi2;
 // fit option
   std::string _peakFitOption; // option flag for root fit
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
    _truncateADC(pset.get<bool>("TruncateADC",true)), 
    _floatPedestal(pset.get<bool>("FloatPedestal",false)), 
    _floatWidth(pset.get<bool>("FloatWidth",false)), 
    _earlyPeak(pset.get<bool>("EarlyPeak",false)),
    _latePeak(pset.get<bool>("LatePeak",false)),
    _findPeaks(pset.get<bool>("FindPeaks",false)),
    _printLevel(pset.get<int>("printLevel",0)),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _debugLevel(pset.get<int>("debugLevel",0)),
    _peakFitOption(pset.get<std::string>("PeakFitOption","QNS")),
    _strawDigis(pset.get<string>("StrawDigis","makeSD"))
  {
    produces<StrawHitCollection>();
    produces<PtrStepPointMCVectorCollection>("StrawHitMCPtr");
    produces<StrawDigiMCCollection>("StrawHitMC");
    if(_printLevel > 0) cout << "In StrawHitsFromStrawDigis constructor " << endl;
  }

  void StrawHitsFromStrawDigis::beginJob(){
    if(_diagLevel > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _shdiag =tfs->make<TTree>("swdiag","StrawHit diagnostics");
      _shdiag->Branch("device",&_shdevice,"device/I");
      _shdiag->Branch("sector",&_shsector,"sector/I");
      _shdiag->Branch("layer",&_shlayer,"layer/I");
      _shdiag->Branch("straw",&_shstraw,"straw/I");
      _shdiag->Branch("charge",&_shq,"charge/I");
      _shdiag->Branch("qstat",&_shqstat,"qstat/I");
      _shdiag->Branch("qchi2",&_shqchi2,"qchi2/F");
    }
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

    // create the peak fit
//    TrkChargeReco::PeakFit pfit(*_strawele);
    TrkChargeReco::FitConfig myconfig;
    myconfig._debug = _debugLevel;
    if(_floatWidth)myconfig.setOption(TrkChargeReco::FitConfig::floatWidth);
    if(_floatPedestal)myconfig.setOption(TrkChargeReco::FitConfig::floatPedestal);
    if(_truncateADC)myconfig.setOption(TrkChargeReco::FitConfig::truncateADC);
    if(_earlyPeak)myconfig.setOption(TrkChargeReco::FitConfig::earlyPeak);
    if(_latePeak)myconfig.setOption(TrkChargeReco::FitConfig::latePeak);
    if(_findPeaks)myconfig.setOption(TrkChargeReco::FitConfig::findPeaks);
    TrkChargeReco::PeakFitRoot pfit(*_strawele,myconfig,_peakFitOption);

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
// fit the ADC waveform to get the charge integral
	StrawDigi::ADCWaveform const& adc = digi.adcWaveform();
	// note: pedestal is being subtracting inside strawele, in the real experiment we will need
	// per-channel version of this FIXME!!!
	TrkChargeReco::PeakFitParams params;
	pfit.process(adc,params);
	if(_debugLevel > 0){
	  cout << "Fit status = " << params._status << " NDF = " << params._ndf << " chisquared " << params._chi2
	  << " Fit charge = " << params._charge << " Fit time = " << params._time << endl;

	}
	// use time division to correct for attenuation FIXME!!
	// the gain should come from a straw-dependent database FIXME!!
	double energy = _strawphys->ionizationEnergy(params._charge/_strawphys->strawGain(2.0,0.0));
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

