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
#include "TrkChargeReco/inc/ComboPeakFitRoot.hh"
#include "TrackerMC/inc/SHInfo.hh"
//CLHEP
#include "CLHEP/Vector/ThreeVector.h"
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
using namespace CLHEP;
namespace mu2e {

  class StrawHitsFromStrawDigis : public art::EDProducer {

  public:
    explicit StrawHitsFromStrawDigis(fhicl::ParameterSet const& pset);
    virtual ~StrawHitsFromStrawDigis(); 
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
    bool _sumADC;
    bool _truncateADC; // model ADC truncation
    bool _floatPedestal; // float pedestal in fit
    bool _floatWidth; // _float width in fit
    bool _earlyPeak, _latePeak, _findPeaks; // additional peak finding
// Diagnostics level.
    int _printLevel, _diagLevel, _debugLevel;
 // fit option
    string _peakFitOption; // option flag for root fit
    unsigned _maxFitIter; //
    // Name of the StrawDigi collection
    string _strawDigis;
    ConditionsHandle<StrawElectronics> _strawele; // models of straw response to stimuli
    ConditionsHandle<StrawPhysics> _strawphys; // models of straw response to stimuli
    TrkChargeReco::PeakFit *_pfit; // peak fitter
// Diagnostics
    TTree* _shdiag;
    SHID _shid; // strawhit ID
    TrkChargeReco::PeakFitParams _peakfit; // result from peak fit
    Float_t _edep, _time, _dt;
    SHMCInfo _shmcinfo; // mc truth info (if available)
    void fillDiagMC(Straw const& straw, StrawDigiMC const& mcdigi);
  };

  StrawHitsFromStrawDigis::StrawHitsFromStrawDigis(fhicl::ParameterSet const& pset) :
    _nbase(pset.get<unsigned>("NumADCBaseline",1)),
    _mbbuffer(pset.get<double>("TimeBuffer",100.0)), // nsec
    _maxdt(pset.get<double>("MaxTimeDifference",8.0)), // nsec
    _singledigi(pset.get<bool>("UseSingleDigis",false)), // use or not single-end digitizations
    _sumADC(pset.get<bool>("SumADC",false)), 
    _truncateADC(pset.get<bool>("TruncateADC",true)), 
    _floatPedestal(pset.get<bool>("FloatPedestal",true)), 
    _floatWidth(pset.get<bool>("FloatWidth",true)), 
    _earlyPeak(pset.get<bool>("EarlyPeak",false)),
    _latePeak(pset.get<bool>("LatePeak",false)),
    _findPeaks(pset.get<bool>("FindPeaks",true)),
    _printLevel(pset.get<int>("printLevel",0)),
    _diagLevel(pset.get<int>("diagLevel",0)),
    _debugLevel(pset.get<int>("debugLevel",0)),
    _peakFitOption(pset.get<string>("PeakFitOption","QNSEX0B")),
    _maxFitIter(pset.get<unsigned>("MaxFitIterations",1)),
    _strawDigis(pset.get<string>("StrawDigis","makeSD")),
    _pfit(0)
  {
    produces<StrawHitCollection>();
    produces<PtrStepPointMCVectorCollection>("StrawHitMCPtr");
    produces<StrawDigiMCCollection>("StrawHitMC");
    if(_printLevel > 0) cout << "In StrawHitsFromStrawDigis constructor " << endl;
  }

  StrawHitsFromStrawDigis::~StrawHitsFromStrawDigis() { delete _pfit; }

  void StrawHitsFromStrawDigis::beginJob(){
    if(_diagLevel > 1){
      art::ServiceHandle<art::TFileService> tfs;
      _shdiag =tfs->make<TTree>("shdiag","StrawHit diagnostics");
      _shdiag->Branch("shid",&_shid); // strawhit ID
      _shdiag->Branch("peakfit",&_peakfit); // ADC waveform fit information
      _shdiag->Branch("edep",&_edep,"edep/F"); // reconstructed deposited energy
      _shdiag->Branch("time",&_time,"time/F"); // reconstructed deposited energy
      _shdiag->Branch("dt",&_dt,"dt/F"); // reconstructed deposited energy
      _shdiag->Branch("mcinfo",&_shmcinfo); // MC truth info
    }
  }

  void StrawHitsFromStrawDigis::beginRun( art::Run& run ){
// create and configure the ADC waveform charge extraction fit
    TrkChargeReco::FitConfig myconfig;
    myconfig._debug = _debugLevel;
    myconfig._maxnit = _maxFitIter;
    if(_floatWidth)myconfig.setOption(TrkChargeReco::FitConfig::floatWidth);
    if(_floatPedestal)myconfig.setOption(TrkChargeReco::FitConfig::floatPedestal);
    if(_truncateADC)myconfig.setOption(TrkChargeReco::FitConfig::truncateADC);
    if(_earlyPeak)myconfig.setOption(TrkChargeReco::FitConfig::earlyPeak);
    if(_latePeak)myconfig.setOption(TrkChargeReco::FitConfig::latePeak);
    if(_findPeaks)myconfig.setOption(TrkChargeReco::FitConfig::findPeaks);
    _strawele = ConditionsHandle<StrawElectronics>("ignored");
    if(_sumADC)
      _pfit = new TrkChargeReco::PeakFit(*_strawele);
    else {
      if(_findPeaks)
	_pfit = new TrkChargeReco::ComboPeakFitRoot(*_strawele,myconfig,_peakFitOption);
      else
	_pfit = new TrkChargeReco::PeakFitRoot(*_strawele,myconfig,_peakFitOption);
    }

    if(_printLevel > 0) cout << "In StrawHitsFromStrawDigis beginRun " << endl;
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
    const Tracker& tracker = getTrackerOrThrow();

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
// fit the ADC waveform to get the charge
	StrawDigi::ADCWaveform const& adc = digi.adcWaveform();
	// note: pedestal is being subtracting inside strawele, in the real experiment we will need
	// per-channel version of this FIXME!!!
	TrkChargeReco::PeakFitParams params;
	_pfit->process(adc,params);
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
// diagnostics
	if(_diagLevel > 1){
	  const Straw& straw = tracker.getStraw( newhit.strawIndex() );
	  _shid = SHID( straw.id() );
	  _peakfit = params;
	  _edep = newhit.energyDep();
	  _time = newhit.time();
	  _dt = newhit.dt();
	  if(mcdigis != 0) fillDiagMC( straw, (*mcdigis)[isd]);
	  _shdiag->Fill();
	}
      }
    }
// put objects into event
    event.put(move(strawHits));
    if(mcptrdigis != 0)event.put(move(mcptrHits),"StrawHitMCPtr");
    if(mchits != 0)event.put(move(mchits),"StrawHitMC");
  }

  void StrawHitsFromStrawDigis::fillDiagMC(Straw const& straw,
    StrawDigiMC const& mcdigi) {
    _shmcinfo._energy = mcdigi.energySum();
    _shmcinfo._trigenergy = mcdigi.triggerEnergySum();
// sum the energy from the explicit trigger particle, and find it's releationship 
    StrawDigi::TDCChannel itdc = StrawDigi::zero;
    if(!mcdigi.hasTDC(itdc)) itdc = StrawDigi::one;
    art::Ptr<StepPointMC> const& spmcp = mcdigi.stepPointMC(itdc);
    art::Ptr<SimParticle> const& spp = spmcp->simParticle();
    _shmcinfo._xtalk = straw.index() != spmcp->strawIndex();
    _shmcinfo._threshenergy = 0.0;
    std::set<art::Ptr<SimParticle> > spptrs; // set of unique particles contributing to this hit
    for(auto imcs = mcdigi.stepPointMCs().begin(); imcs!= mcdigi.stepPointMCs().end(); ++ imcs){
    // if the simParticle for this step is the same as the one which fired the discrim, add the energy
      if( (*imcs)->simParticle() == spp )
	_shmcinfo._threshenergy += (*imcs)->eDep();
      spptrs.insert((*imcs)->simParticle()); 
    }
    _shmcinfo._nmcpart = spptrs.size();
    _shmcinfo._pdg = spp->pdgId();
    _shmcinfo._proc = spp->creationCode();
    if(spp->genParticle().isNonnull())
      _shmcinfo._gen = spp->genParticle()->generatorId().id();
    else
      _shmcinfo._gen=-1;
    Hep3Vector mcsep = spmcp->position()-straw.getMidPoint();
    Hep3Vector dir = spmcp->momentum().unit();
    _shmcinfo._mom = spmcp->momentum().mag();
    Hep3Vector mcperp = (dir.cross(straw.getDirection())).unit();
    double dperp = mcperp.dot(mcsep);
    _shmcinfo._dperp = fabs(dperp);
    _shmcinfo._ambig = dperp > 0 ? -1 : 1; // follow TrkPoca convention
    _shmcinfo._len = mcsep.dot(straw.getDirection());

  }
}
using mu2e::StrawHitsFromStrawDigis;
DEFINE_ART_MODULE(StrawHitsFromStrawDigis);

