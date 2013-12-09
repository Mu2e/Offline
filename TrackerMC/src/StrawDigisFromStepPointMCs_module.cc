//
// This module transforms StepPointMC objects into StrawDigi objects
// It also builds the truth match map
//
// $Id: StrawDigisFromStepPointMCs_module.cc,v 1.3 2013/12/09 05:09:50 brownd Exp $
// $Author: brownd $ 
// $Date: 2013/12/09 05:09:50 $
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
#include "SeedService/inc/SeedService.hh"
#include "cetlib/exception.h"
// conditions
#include "ConditionsService/inc/ConditionsHandle.hh"
#include "ConditionsService/inc/TrackerCalibrations.hh"
#include "ConditionsService/inc/AcceleratorParams.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TTrackerGeom/inc/TTracker.hh"
#include "ConfigTools/inc/ConfigFileLookupPolicy.hh"
#include "HitMakers/inc/DeadStrawList.hh"
// data
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// MC structures
#include "TrackerMC/inc/StrawHitletSequencePair.hh"
#include "TrackerMC/inc/StrawWaveform.hh"
//CLHEP
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"
// root 
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TTree.h"
#// C++
#include <map>
using namespace std;

namespace mu2e {

  struct IonCluster {  // ion charge cluster before drift or amplification
    CLHEP::Hep3Vector _pos; // position of this cluster
    double _charge; // charge of this cluster, in pC.  Note: this is pre-gain!!!
    IonCluster(CLHEP::Hep3Vector const& pos, double charge): _pos(pos),_charge(charge) {} 
  };

  struct WireCharge { // charge at the wire after drift
    double _charge; // charge at the wire, in units of pC
    double _time; // relative time at the wire, relative to ionzation time (ns)
    double _wpos; // position long the wire, WRT the wire center, signed by the wire direction
  };

  struct WireEndCharge { // charge at one end of the wire after propagation
    double _charge; // charge at the wire, in units of pC
    double _time; // time at the wire end, relative to the time the charge reached the wire (ns)
    double _wdist; // propagation distance from the point of collection to the end
  };

  class StrawDigisFromStepPointMCs : public art::EDProducer {

  public:

    typedef map<StrawIndex,StrawHitletSequencePair> StrawHitletMap;  // hitlets by straw 
    typedef vector<art::Ptr<StepPointMC> > StrawSPMCPV; // vector of associated StepPointMCs for a single straw/particle
    typedef list<WFX> WFXList;
    typedef vector<WFXList::const_iterator> WFXP;

    explicit StrawDigisFromStepPointMCs(fhicl::ParameterSet const& pset);
    // Accept compiler written d'tor.

    virtual void beginJob();
    virtual void beginRun( art::Run& run );
    void produce( art::Event& e);

  private:

    // Diagnostics level.
    int _diagLevel;

    // Limit on number of events for which there will be full printout.
    int _maxFullPrint;

    // Name of the tracker StepPoint collection
    string _trackerStepPoints;

    // Parameters
    bool   _addXtalk; // should we add cross talk hits?
    bool   _addNoise; // should we add noise hits?
    string _g4ModuleLabel;  // Name of the module that made these hits.
    double _minsteplen; // minimum step size for ionization cluster
// the following parmeters describe straw drift and should be in a separate calibration object FIXME!!
    double _EIonize; // Geant energy of each ionization (MeV)
    double _QIonize; // charge of a single ionization (=e, pC)
    double _gasgain; // avalanche gain
    double _attlen; // attenuation length of charge down the wire
    double _vprop; // propagation time of signals down the wire, mm/nsec
    double _vdrift; // electron drift velocity
    double _vdrifterr; // drift velocity error
    double _mbtime; // period of 1 microbunch
    double _mbbuffer; // buffer on that for ghost hits (wrapping)
    double _vthresh; // threshold voltage for electronics discriminator
    StrawElectronics _strawele; // models of straw response to stimuli
    // Random number distributions
    CLHEP::RandGaussQ _gaussian;
    CLHEP::RandFlat _randflat;
    // A category for the error logger.
    const string _messageCategory;
    // Give some informationation messages only on the first event.
    bool _firstEvent;
    // List of dead straws as a parameter set; needed at beginRun time.
    fhicl::ParameterSet _deadStraws;
    // Access the live/dead status of the straw.
    DeadStrawList _strawStatus;
// diagnostics
    TTree* _sddiag;
    Int_t _device, _sector, _layer, _straw;
    Float_t _time0, _time1, _mctime;
    vector<unsigned> _adc;
    Int_t _t0, _t1;
    TH1F* _waveforms;
//  helper functions
    void fillHitletMap(art::Event const& event, StrawHitletMap & hmap);
    void addStep(art::Ptr<StepPointMC>& spmcptr, Straw const& straw, StrawHitletSequencePair& shsp);
    void divideStep(StepPointMC const& step, vector<IonCluster>& clusters);
    void driftCluster(Straw const& straw, IonCluster const& cluster, WireCharge& wireq);
    void propagateCharge(Straw const& straw, WireCharge const& wireq, StrawEnd end, WireEndCharge& weq);
    double microbunchTime(double globaltime) const;
    void addGhosts(StrawHitlet const& hitlet,StrawHitletSequence& shs); 
    void addCrosstalk(StrawHitletMap& hmap);
    void addNoise(StrawHitletMap& hmap);
    void findThresholdCrossings(StrawHitletSequencePair const& hsp, WFXList& xings);
    void fillDigis(WFXList const& xings,StrawWaveform const& wf, StrawDigiCollection* digis);
    void createDigi(WFXP const& xpair, StrawWaveform const& primarywf, StrawDigiCollection* digis);
 
    // the following should be delegated to a straw physics description object, FIXME!!
    double strawGain(Straw const& straw, double driftdist, double driftphi) const;
    void distanceToTime(Straw const& straw, double driftdist, double driftphi,
	double& drifttime, double& drifttimeerr) const;
    StrawEnd primaryEnd(StrawIndex straw_id) const;
  };

  StrawDigisFromStepPointMCs::StrawDigisFromStepPointMCs(fhicl::ParameterSet const& pset) :

    // Parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _maxFullPrint(pset.get<int>("maxFullPrint",5)),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _addXtalk(pset.get<bool>("addCrossTalkHits",false)),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _minsteplen(pset.get<double>("MinimumIonClusterStep",0.5)),
    _EIonize(pset.get<double>("EnergyPerIonization",27.0e-6)), // 27 ev
    _QIonize(pset.get<double>("ChargePerIonization",1.6e-7)), // e, pC
    _gasgain(pset.get<double>("GasGain",1.0e5)),
    _attlen(pset.get<double>("PropagationAttentuationLength",1000.0)), // mm
    _vdrift(pset.get<double>("DriftVelocity",0.05)), // mm/nsec
    _vdrifterr(pset.get<double>("DriftVelocity",0.001)), // mm/nsec
    _mbbuffer(pset.get<double>("TimeFoldingBuffer",200.0)), // nsec
    _vthresh(pset.get<double>("DiscriminatorThreshold",50.0)), //mVolt
    _strawele(pset.get<fhicl::ParameterSet>("StrawElectronics",fhicl::ParameterSet())),

    // Random number distributions
    _gaussian( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),
    _randflat( createEngine( art::ServiceHandle<SeedService>()->getSeed() ) ),

    _messageCategory("HITS"),

    // Control some information messages.
    _firstEvent(true),

    _deadStraws(pset.get<fhicl::ParameterSet>("deadStrawList", fhicl::ParameterSet())),
    _strawStatus(_deadStraws) 
    {
 
// Tell the framework what we make.
      produces<StrawDigiCollection>();
// should also produce a MC truth map from StrawDigis to StepPointMCs:  FIXME!!
    }

  void StrawDigisFromStepPointMCs::beginJob(){
    if(_diagLevel > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _sddiag =tfs->make<TTree>("sddiag","StrawDigi diagnostics");
      _sddiag->Branch("device",&_device,"device/I");
      _sddiag->Branch("sector",&_sector,"sector/I");
      _sddiag->Branch("layer",&_layer,"layer/I");
      _sddiag->Branch("straw",&_straw,"straw/I");
      _sddiag->Branch("time0",&_time0,"time0/F");
      _sddiag->Branch("time1",&_time1,"time1/F");
      _sddiag->Branch("t0",&_t0,"t0/I");
      _sddiag->Branch("t1",&_t1,"t1/I");
      _sddiag->Branch("adc",&_adc);
      _sddiag->Branch("mctime",&_mctime,"mctime/F");
    }
    if(_diagLevel > 1){
// record a few waveforms

    }
  }

  void StrawDigisFromStepPointMCs::beginRun( art::Run& run ){
    _strawStatus.reset(_deadStraws);
  }

  void
  StrawDigisFromStepPointMCs::produce(art::Event& event) {
    if ( _diagLevel > 1 ) cout << "StrawDigisFromStepPointMCs: produce() begin; event " << event.id().event() << endl;
    static int ncalls(0);
    ++ncalls;
// update conditions caches.  The conditions handles themselves should be data members FIXME!!
    ConditionsHandle<AcceleratorParams> accPar("ignored");
    _mbtime = accPar->deBuncherPeriod;
   // Containers to hold the output information.
    unique_ptr<StrawDigiCollection> digis(new StrawDigiCollection);
    // Handle to the conditions service
    ConditionsHandle<TrackerCalibrations> trackerCalibrations("ignored");
// create the StrawHitlet map
    StrawHitletMap hmap;
// fill this from the event
    fillHitletMap(event,hmap);
// introduce cross-talk hitlets
    if(_addXtalk)addCrosstalk(hmap);
// add noise hitlets
    if(_addNoise)addNoise(hmap);
// loop over the hitlet sequences
    for(auto ihsp=hmap.begin();ihsp!= hmap.end();++ihsp){
// find the threshold crossing points along this hitlet sequence
      WFXList xings;
      findThresholdCrossings(ihsp->second,xings);
// convert the crossing points into digis, and add them to the event data
      if(xings.size() > 0){
// instantiate a waveform for the primary end of this straw
	StrawEnd primaryend = primaryEnd(ihsp->second.strawIndex());
	StrawWaveform primarywf(ihsp->second.hitletSequence(primaryend),_strawele);
// create digis
	fillDigis(xings,primarywf,digis.get());
      }
    }
// store the digis in the event
    event.put(move(digis));
    if ( _diagLevel > 1 ) cout << "StrawDigisFromStepPointMCs: produce() end" << endl;
    // Done with the first event; disable some messages.
    _firstEvent = false;
  } // end produce

  void
  StrawDigisFromStepPointMCs::fillHitletMap(art::Event const& event, StrawHitletMap & hmap){
    // get conditions
    const Tracker& tracker = getTrackerOrThrow();
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    // Get all of the tracker StepPointMC collections from the event:
    typedef std::vector< art::Handle<StepPointMCCollection> > HandleVector;
    // This selector will select only data products with the given instance name.
    art::ProductInstanceNameSelector selector("tracker");
    HandleVector stepsHandles;
    event.getMany( selector, stepsHandles);
    // Informational message on the first event.
    if ( _firstEvent ) {
      mf::LogInfo log(_messageCategory);
      log << "StrawDigisFromStepPointMCs::fillHitMap will use StepPointMCs from: \n";
      for ( HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end();
	  i != e; ++i ){
	art::Provenance const& prov(*(i->provenance()));
	log  << "   " << prov.branchName() << "\n";
      }
    }
    // Loop over StepPointMC collections
    for ( HandleVector::const_iterator ispmcc=stepsHandles.begin(), espmcc=stepsHandles.end();ispmcc != espmcc; ++ispmcc ){
      art::Handle<StepPointMCCollection> const& handle(*ispmcc);
      StepPointMCCollection const& steps(*handle);
 // Loop over the StepPointMCs in this collection
      for (size_t ispmc =0; ispmc<steps.size();++ispmc){
      // find straw index
        StrawIndex const & straw_id = steps[ispmc].strawIndex();
	// lookup straw here, to avoid having to find the tracker for every step
        Straw const& straw = tracker.getStraw(straw_id);
	// Skip dead straws.
	if ( _strawStatus.isDead(straw_id)) continue;
	// Skip steps that occur in the deadened region near the end of each wire.
	double wpos = fabs((steps[ispmc].position()-straw.getMidPoint()).dot(straw.getDirection()));
	if(wpos >  straw.getDetail().activeHalfLength())continue;
	// create ptr to MC truth, used for references
	art::Ptr<StepPointMC> spmcptr(handle,ispmc);
	// create a hitlet from this step, and add it to the hitlet map
	addStep(spmcptr,straw,hmap[straw_id]);
      }
    }
  }

  void
  StrawDigisFromStepPointMCs::addStep(art::Ptr<StepPointMC>& spmcptr,
      Straw const& straw,  
      StrawHitletSequencePair& shsp) {
    StepPointMC const& step = *spmcptr;
    StrawIndex const & straw_id = step.strawIndex();
    // Subdivide the StepPointMC into ionization clusters
    vector<IonCluster> clusters;
    divideStep(step,clusters);
    // drift these clusters to the wire, and record the charge at the wire
    for(auto iclu=clusters.begin(); iclu != clusters.end(); ++iclu){
      WireCharge wireq;
      driftCluster(straw,*iclu,wireq);
// propagate this charge to each end of the wire
      for(size_t iend=0;iend<2;++iend){
	StrawEnd end(static_cast<StrawEnd::strawend>(iend));
	// compute the longitudinal propagation effects
	WireEndCharge weq;
	propagateCharge(straw,wireq,end,weq);
	// compute the total time, modulo the microbunch
	double gtime = step.time()+weq._time + weq._time;
	double htime = microbunchTime(gtime);
	// create the hitlet
	StrawHitlet hitlet(StrawHitlet::primary,straw_id,end,htime,weq._charge,weq._wdist,spmcptr);
	// add the hitlets to the appropriate sequence.
	shsp.hitletSequence(end).insert(hitlet);
	// if required, add a 'ghost' copy of this hitlet
	addGhosts(hitlet,shsp.hitletSequence(end));
      }
    }
  }

  void
  StrawDigisFromStepPointMCs::divideStep(StepPointMC const& step,
      vector<IonCluster>& clusters) {
// if the step is already small enough, don't subdivide
    if(step.stepLength() < _minsteplen) {
      IonCluster cluster(step.position(),_QIonize*step.ionizingEdep()/_EIonize); 
      clusters.push_back(cluster);
    } else {
// subdivide into units of fundamental charge, but no smaller than the smallest step size
      unsigned maxndiv = static_cast<unsigned>(ceil(step.stepLength()/_minsteplen));
      static unsigned one(1);
      unsigned ndiv = min(maxndiv,max(one,static_cast<unsigned>(rint(step.ionizingEdep()/_EIonize))));
// generate random points for each ionization
      vector<double> lengths(ndiv,0.0);
// the following assumes StepPointMC::position is at the begining of the step
      _randflat.shootArray(ndiv,lengths.data(),0.0,step.stepLength());
// sort these.  I'm not sure this is necessary CHECKME!!!
      sort(lengths.begin(),lengths.end());
// make clusters for each
      CLHEP::Hep3Vector dir = step.momentum().unit();
      double qdiv = _QIonize*step.ionizingEdep()/(_EIonize*ndiv); 
      for(auto ilen=lengths.begin();ilen!=lengths.end();++ilen){
	CLHEP::Hep3Vector pos = step.position() + (*ilen)*dir;
	IonCluster cluster(pos,qdiv);
	clusters.push_back(cluster);
      }
    }
  }

  void
  StrawDigisFromStepPointMCs::driftCluster(Straw const& straw,
      IonCluster const& cluster, WireCharge& wireq) {
// Compute the vector from the cluster to the wire
    CLHEP::Hep3Vector cpos = cluster._pos-straw.getMidPoint();
    // drift distance perp to wire, and angle WRT magnetic field (for Lorentz effect)
    double dd = cpos.perp(straw.getDirection());
    double dphi = cpos.azimAngle(straw.getDirection());
    double gain = strawGain(straw,dd,dphi);
    // smear the charge by the gas gain statistics 
    double dgain = _gaussian.shoot(0.0,sqrt(gain));
    wireq._charge = cluster._charge*(gain+dgain);
    // time is smeared proportional to drift time
    double dtime, dtimeerr;
    distanceToTime(straw,dd,dphi,dtime,dtimeerr);
    wireq._time = _gaussian.shoot(dtime,dtimeerr);
    // position along wire
    wireq._wpos = cpos.dot(straw.getDirection());
  }

  double
  StrawDigisFromStepPointMCs::strawGain(Straw const& straw, double driftdist, double driftphi) const {
      // someday need a physics model and calibration constants here, FIXME!!!
      return _gasgain;
  }

  void
  StrawDigisFromStepPointMCs::distanceToTime(Straw const& straw, double driftdist, double driftphi,
    double& drifttime, double& drifttimeerr) const {
// again need a better model here, FIXME!!
  // I'm not sure how expensive this call is: it could be performed once/event then passed in FIXME!!
//    ConditionsHandle<TrackerCalibrations> tcal("ignored");
//    static T2D t2d;
//    static CLHEP::Hep3Vector tdir(0.0,0.0,1.0);
//    tcal->TimeToDistance(straw().index(),tdrift,tdir,t2d);
    // I'm using t2d to compute dfit time, which isn't right, since t2d refers to tracks and this is for clusters
    // FIXME!!
//    drifttime =  driftdist/t2d._vdrift;
    drifttime =  driftdist/_vdrift;
    drifttimeerr = driftdist/_vdrifterr;
  }

  void
  StrawDigisFromStepPointMCs::propagateCharge(Straw const& straw,
      WireCharge const& wireq, StrawEnd end, WireEndCharge& weq) {
  // I'm not sure how expensive this call is: it could be performed once/event then cached or passed in FIXME!!
    ConditionsHandle<TrackerCalibrations> tcal("ignored");
    double vwire = tcal->SignalVelocity(straw.index());
    // compute distance to the appropriate end
    double wlen = straw.getDetail().halfLength(); // use the full length, not the active length
  // NB: the following logic assumes the straw direction points in increasing azimuth.  FIXME!
    if(end == StrawEnd::plus)
      weq._wdist = wlen - wireq._wpos;
    else
      weq._wdist = wlen + wireq._wpos;
    // split the charge, and attenuate it according to the distance
    weq._charge = 0.5*wireq._charge*exp(-weq._wdist/_attlen);
    // linear time propagation.  Dispersion is handled elsewhere
    weq._time = weq._wdist/vwire;
  }

  double
  StrawDigisFromStepPointMCs::microbunchTime(double globaltime) const {
    return fmod(globaltime,_mbtime);
  }

  void
  StrawDigisFromStepPointMCs::addGhosts(StrawHitlet const& hitlet,StrawHitletSequence& shs) {
    if(hitlet.time() < _mbbuffer)
      shs.insert(StrawHitlet(hitlet,_mbtime));
    if(_mbtime-hitlet.time() < _mbbuffer)
      shs.insert(StrawHitlet(hitlet,-_mbtime));
  }

  void
  StrawDigisFromStepPointMCs::findThresholdCrossings(StrawHitletSequencePair const& hsp, WFXList& xings){
    xings.clear();
    // loop over the ends of this straw
    for(size_t iend=0;iend<2;++iend){
      StrawEnd end(static_cast<StrawEnd::strawend>(iend));
    // convert the hitlet list from this end to waveforms.  This is a light-weight process
      StrawWaveform swf(hsp.hitletSequence(end),_strawele);
// iterate sequentially over hitlets inside the sequence
      WFX wfx(swf); // start at the begining of the microbunch, at the begining of the wavelets
      // find where the waveform xings threshold.  This updates the waveform sample
      // Skip any points outside the microbunch readout window
      while(swf.crossesThreshold(_vthresh,wfx) && wfx._time < _mbtime){
      // keep these in time-order
	auto iwfxl = xings.begin();
	while(iwfxl != xings.end() && iwfxl->_time < wfx._time)++iwfxl;
	xings.insert(iwfxl,wfx);
      }
    }
  }

  void
  StrawDigisFromStepPointMCs::fillDigis(WFXList const& xings, StrawWaveform const& primarywf,
    StrawDigiCollection* digis){
// loop over crossings
    auto iwfxl = xings.begin();
    while(iwfxl!= xings.end()){
      WFXP xpair(1,iwfxl);
 // associate adjacent crossing if they are on opposite ends within the maximum propagation time difference
      auto jwfxl = ++iwfxl;
      if(jwfxl != xings.end() &&
	  iwfxl->_ihitlet->strawEnd() == jwfxl->_ihitlet->strawEnd() &&
	  jwfxl != xings.end() && _strawele.combineEnds(iwfxl->_time,jwfxl->_time)) {
	xpair.push_back(jwfxl);
	iwfxl = ++jwfxl;
      } else
	iwfxl = jwfxl;
// create a digi from this pair or singleton
      createDigi(xpair,primarywf,digis);
      if(_diagLevel > 0){
	const Tracker& tracker = getTrackerOrThrow();
	const Straw& straw = tracker.getStraw( primarywf.hitlets().strawIndex() );
	_device = straw.id().getDevice();
	_sector = straw.id().getSector();
	_layer = straw.id().getLayer();
	_straw = straw.id().getStraw();
	for(auto ixp=xpair.begin();ixp!=xpair.end();++ixp){
	  if((*ixp)->_ihitlet->strawEnd() == primarywf.hitlets().strawEnd())
	    _time0 =(*ixp)->_time;
	  else
	    _time1 =(*ixp)->_time;
	}
	_mctime = xpair[0]->_ihitlet->stepPointMC()->time();
	StrawDigi const& digi = digis->back();
	_t0 = digi.TDC(StrawDigi::zero);
	_t1 = digi.TDC(StrawDigi::one);
	_adc.clear();
	for(auto iadc=digi.adcWaveform().begin();iadc!=digi.adcWaveform().end();++iadc){
	  _adc.push_back(*iadc);
	}
      }
    }
  }

  void StrawDigisFromStepPointMCs::createDigi(WFXP const& xpair, StrawWaveform const& primarywf,
    StrawDigiCollection* digis){
// storage for MC match can be more than 1 StepPointMCs
    set<art::Ptr<StepPointMC>> mcmatch;
// initialize the float variables that we later digitize
    array<double,2> xtimes = {0.0,0.0};
    vector<double> wf(_strawele.nADCSamples(),0.0);
// loop over the associated crossings
    for(auto iwfx = xpair.begin();iwfx!= xpair.end();++iwfx){
      WFX const& wfx = **iwfx;
    // primary end times is index 0, other end is index 1
      size_t index(1);
      if(wfx._ihitlet->strawEnd() == primarywf.hitlets().strawEnd()){
	index = 0;
// sample the waveform at the primary end
	vector<double> adctimes;
	_strawele.adcTimes(wfx._time,adctimes);
	primarywf.sampleWaveform(adctimes,wf);
      }
      xtimes[index] = wfx._time; // non-primary end is TDC 1;
  // record MC match if it isn't already recorded
      mcmatch.insert(wfx._ihitlet->stepPointMC());
    }
// digitize
    StrawDigi::ADCWaveform adc;
    _strawele.digitizeWaveform(wf,adc);
    StrawDigi::TDCValues tdc;
    _strawele.digitizeTimes(xtimes,tdc);
// create the digi from this
    digis->push_back(StrawDigi(primarywf.hitlets().strawIndex(),tdc,adc));
// fill map entry to record association of this digi with StepPointMC.
// FIXME!!!!
  }

  StrawEnd
  StrawDigisFromStepPointMCs::primaryEnd(StrawIndex straw_id) const {
    // Here, Im assuming the plus end is the primary, whereas in the real detector the primary
    // ends will alternate and needs to be looked up in an electronics map FIXME!!!!!
    return StrawEnd(StrawEnd::plus);
  }

  // functions that need implementing:: FIXME!!!!!!
  void  StrawDigisFromStepPointMCs::addCrosstalk(StrawHitletMap& hmap) {
    // loop over straws
  // find nearby straws with non-zero coupling constants
  // create secondary hitlets from the primary ones, and insert them into the coulpled straw hitlet sequence
  }
  void StrawDigisFromStepPointMCs::addNoise(StrawHitletMap& hmap){
  // create random noise hitlets and add them to the sequences of random straws.
  }

} // end namespace mu2e

using mu2e::StrawDigisFromStepPointMCs;
DEFINE_ART_MODULE(StrawDigisFromStepPointMCs)


