//
// This module transforms StepPointMC objects into StrawDigi objects
// It also builds the truth match map
//
// $Id: StrawDigisFromStepPointMCs_module.cc,v 1.6 2013/12/12 00:39:24 brownd Exp $
// $Author: brownd $ 
// $Date: 2013/12/12 00:39:24 $
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
// utiliities
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
//#include "HitMakers/inc/DeadStrawList.hh"
// data
#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
// MC structures
#include "TrackerMC/inc/StrawHitletSequencePair.hh"
#include "TrackerMC/inc/StrawElectronics.hh"
#include "TrackerMC/inc/StrawWaveform.hh"
//CLHEP
#include "CLHEP/Random/RandGaussQ.h"
#include "CLHEP/Random/RandFlat.h"
// root 
#include "TMath.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
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
    double _dd;  // transverse distance drifted to the wrie
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
    unsigned _maxhist;
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
    art::RandomNumberGenerator::base_engine_t& _engine;
    CLHEP::RandGaussQ _gaussian;
    CLHEP::RandFlat _randflat;
    // A category for the error logger.
    const string _messageCategory;
    // Give some informationation messages only on the first event.
    bool _firstEvent;
    unsigned _event;
    // List of dead straws as a parameter set; needed at beginRun time.
//    DeadStrawList _strawStatus;
// diagnostics
    TTree* _sdiag;
    Int_t _sdevice, _ssector, _slayer, _sstraw;
    Int_t _nhitlet;
    Int_t _nsteppoint;
    Int_t _npart;
    Float_t _tmin, _tmax;
    TTree* _sddiag;
    Int_t _sddevice, _sdsector, _sdlayer, _sdstraw;
    Int_t _nend, _nstep;
    Float_t _xtime0, _xtime1, _htime0, _htime1, _charge0, _charge1, _ddist0, _ddist1, _wdist0, _wdist1;
    Float_t _vmax;
    Float_t _mctime, _mcenergy, _mcdca;
    vector<unsigned> _adc;
    Int_t _tdc0, _tdc1;
//    vector<TGraph*> _waveforms;
    vector<TH1F*> _waveforms;
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

// diagnostic parameters
    _diagLevel(pset.get<int>("diagLevel",0)),
    _maxhist(pset.get<unsigned>("MaxHist",2)),
    // Parameters
    _maxFullPrint(pset.get<int>("maxFullPrint",5)),
    _trackerStepPoints(pset.get<string>("trackerStepPoints","tracker")),
    _addXtalk(pset.get<bool>("addCrossTalkHits",false)),
    _g4ModuleLabel(pset.get<string>("g4ModuleLabel")),
    _minsteplen(pset.get<double>("MinimumIonClusterStep",0.1)), // 100 microns
    _EIonize(pset.get<double>("EnergyPerIonization",27.0e-6)), // 27 ev
    _QIonize(pset.get<double>("ChargePerIonization",1.6e-7)), // e, pC
    _gasgain(pset.get<double>("GasGain",3.0e4)),
    _attlen(pset.get<double>("PropagationAttentuationLength",10000.0)), // mm
    _vdrift(pset.get<double>("DriftVelocity",0.05)), // mm/nsec
    _vdrifterr(pset.get<double>("DriftVelocity",0.001)), // mm/nsec
    _mbbuffer(pset.get<double>("TimeFoldingBuffer",200.0)), // nsec
    _vthresh(pset.get<double>("DiscriminatorThreshold",5.0)), //mVolt
    _strawele(pset.get<fhicl::ParameterSet>("StrawElectronics",fhicl::ParameterSet())),
    // Random number distributions
    _engine(createEngine( art::ServiceHandle<SeedService>()->getSeed())),
    _gaussian( _engine ),
    _randflat( _engine ),

    _messageCategory("HITS"),

    // Control some information messages.
    _firstEvent(true)
//    _strawStatus(pset.get<fhicl::ParameterSet>("deadStrawList", fhicl::ParameterSet()))
    {
// Tell the framework what we make.
      produces<StrawDigiCollection>();
// should also produce a MC truth map from StrawDigis to StepPointMCs:  FIXME!!
    }

  void StrawDigisFromStepPointMCs::beginJob(){
    if(_diagLevel > 0){
      art::ServiceHandle<art::TFileService> tfs;
      _sdiag =tfs->make<TTree>("sdiag","StrawDigi diagnostics");
      _sdiag->Branch("device",&_sdevice,"device/I");
      _sdiag->Branch("sector",&_ssector,"sector/I");
      _sdiag->Branch("layer",&_slayer,"layer/I");
      _sdiag->Branch("straw",&_sstraw,"straw/I");
      _sdiag->Branch("nhitlet",&_nhitlet,"straw/I");
      _sdiag->Branch("nstep",&_nsteppoint,"straw/I");
      _sdiag->Branch("npart",&_npart,"straw/I");
      _sdiag->Branch("tmin",&_tmin,"tmin/F");
      _sdiag->Branch("tmax",&_tmax,"tmax/F");

      _sddiag =tfs->make<TTree>("sddiag","StrawDigi diagnostics");
      _sddiag->Branch("device",&_sddevice,"device/I");
      _sddiag->Branch("sector",&_sdsector,"sector/I");
      _sddiag->Branch("layer",&_sdlayer,"layer/I");
      _sddiag->Branch("straw",&_sdstraw,"straw/I");
      _sddiag->Branch("nend",&_nend,"nend/I");
      _sddiag->Branch("nstep",&_nstep,"nstep/I");
      _sddiag->Branch("xtime0",&_xtime0,"xtime0/F");
      _sddiag->Branch("xtime1",&_xtime1,"xtime1/F");
      _sddiag->Branch("htime0",&_htime0,"htime0/F");
      _sddiag->Branch("htime1",&_htime1,"htime1/F");
      _sddiag->Branch("charge0",&_charge0,"charge0/F");
      _sddiag->Branch("charge1",&_charge1,"charge1/F");
      _sddiag->Branch("wdist0",&_wdist0,"wdist0/F");
      _sddiag->Branch("wdist1",&_wdist1,"wdist1/F");
      _sddiag->Branch("ddist0",&_ddist0,"ddist0/F");
      _sddiag->Branch("ddist1",&_ddist1,"ddist1/F");
      _sddiag->Branch("tdc0",&_tdc0,"tdc0/I");
      _sddiag->Branch("tdc1",&_tdc1,"tdc1/I");
      _sddiag->Branch("adc",&_adc);
      _sddiag->Branch("vmax",&_vmax,"vmax/F");
      _sddiag->Branch("mctime",&_mctime,"mctime/F");
      _sddiag->Branch("mcenergy",&_mcenergy,"mcenergy/F");
      _sddiag->Branch("mcdca",&_mcdca,"mcdca/F");
    }
  }

  void StrawDigisFromStepPointMCs::beginRun( art::Run& run ){
//    _strawStatus.reset(_deadStraws);
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
    if ( _diagLevel > 2 ) cout << "StrawDigisFromStepPointMCs: produce() end" << endl;
    // Done with the first event; disable some messages.
    _firstEvent = false;
    ++_event;
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
//	if ( _strawStatus.isDead(straw_id)) continue;
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
	double gtime = step.time()+wireq._time + weq._time;
	double htime = microbunchTime(gtime);
	// create the hitlet
	StrawHitlet hitlet(StrawHitlet::primary,straw_id,end,htime,weq._charge,wireq._dd,weq._wdist,spmcptr);
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
    // distance to time error is broken!  don't use for now, FIXME!!!
//    wireq._time = _gaussian.shoot(dtime,dtimeerr);
    wireq._time = dtime; 
    wireq._dd = dd;
    // position along wire
    // need to add Lorentz effects, FIXME!!!
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
	while(iwfxl != xings.end() && iwfxl->_time < wfx._time)
	  ++iwfxl;
	xings.insert(iwfxl,wfx);
	// skip to the next hitlet, and insure a minimum time buffer between crossings
	++(wfx._ihitlet);
	wfx._time += _strawele.deadTime();
      }
    }
  }

  void
  StrawDigisFromStepPointMCs::fillDigis(WFXList const& xings, StrawWaveform const& primarywf,
      StrawDigiCollection* digis){
    const Tracker& tracker = getTrackerOrThrow();
    if(_diagLevel > 0){
      const Straw& straw = tracker.getStraw( primarywf.hitlets().strawIndex() );
      _sdevice = straw.id().getDevice();
      _ssector = straw.id().getSector();
      _slayer = straw.id().getLayer();
      _sstraw = straw.id().getStraw();
      HitletList const& hitlets = primarywf.hitlets().hitletList();
      _nhitlet = hitlets.size();
      set<art::Ptr<StepPointMC> > steps;
      set<art::Ptr<SimParticle> > parts;
      for(auto ihitl=hitlets.begin();ihitl!=hitlets.end();++ihitl){
	if(ihitl->stepPointMC().isNonnull()){
	  steps.insert(ihitl->stepPointMC());
	  parts.insert(ihitl->stepPointMC()->simParticle());
	}
	_nsteppoint = steps.size();
	_npart = parts.size();
      }
      _tmin = hitlets.begin()->time();
      _tmax = hitlets.rbegin()->time();
      _sdiag->Fill();
      if(_diagLevel > 2 && _event < _maxhist){
// histogram the waveforms
	const double tstep(0.5); // 0.5
	const double nfall(3.0); // 4 lambda past last fall time
	double tstart = hitlets.begin()->time()-tstep;
	double tfall = _strawele.fallTime();
	double tend = hitlets.rbegin()->time() + nfall*tfall;
	vector<double> times, volts;
	double t = tstart;
	while(t<tend){
	  times.push_back(t);
	  t += tstep;
	}
	primarywf.sampleWaveform(times,volts);
	art::ServiceHandle<art::TFileService> tfs;
	char name[60];
	snprintf(name,60,"SWF%i_%i",primarywf.hitlets().strawIndex().asInt(),_event);
	TH1F* wf = tfs->make<TH1F>(name,name,volts.size(),times.front(),times.back());
	for(size_t ibin=0;ibin<times.size();++ibin)
	  wf->SetBinContent(ibin+1,volts[ibin]);
	_waveforms.push_back(wf);
      }
    }
    // loop over crossings
    auto iwfxl = xings.begin();
    while(iwfxl!= xings.end()){
      WFXP xpair(1,iwfxl);
 // associate adjacent crossing if they are on opposite ends within the maximum propagation time difference
      auto jwfxl = iwfxl; ++jwfxl;
      if(jwfxl != xings.end() &&
	  iwfxl->_ihitlet->strawEnd() != jwfxl->_ihitlet->strawEnd() &&
	  _strawele.combineEnds(iwfxl->_time,jwfxl->_time)) {
	xpair.push_back(jwfxl);
	iwfxl = jwfxl;
      } 
      ++iwfxl;
// create a digi from this pair or singleton
      createDigi(xpair,primarywf,digis);
  // diagnostics
      if(_diagLevel > 1){
	const Straw& straw = tracker.getStraw( primarywf.hitlets().strawIndex() );
	_sddevice = straw.id().getDevice();
	_sdsector = straw.id().getSector();
	_sdlayer = straw.id().getLayer();
	_sdstraw = straw.id().getStraw();
	_xtime0 = _xtime1 = -1000.0;
	_htime0 = _htime1 = -1000.0;
	_charge0 = _charge1 = -1000.0;
	_wdist0 = _wdist1 = -1000.0;
	_ddist0 = _ddist1 = -1000.0;
	_nend = xpair.size();
	for(auto ixp=xpair.begin();ixp!=xpair.end();++ixp){
	  if((*ixp)->_ihitlet->strawEnd() == primarywf.hitlets().strawEnd()) {
	    _xtime0 =(*ixp)->_time;
	    _htime0 = (*ixp)->_ihitlet->time();
	    _charge0 = (*ixp)->_ihitlet->charge();
	    _ddist0 = (*ixp)->_ihitlet->driftDistance();
	    _wdist0 = (*ixp)->_ihitlet->wireDistance();
	  } else {
	    _xtime1 =(*ixp)->_time;
	    _htime1 = (*ixp)->_ihitlet->time();
	    _charge1 = (*ixp)->_ihitlet->charge();
	    _ddist1 = (*ixp)->_ihitlet->driftDistance();
	    _wdist1 = (*ixp)->_ihitlet->wireDistance();
	  }
	}
	if(xpair.size() < 2 || xpair[0]->_ihitlet->stepPointMC() == xpair[1]->_ihitlet->stepPointMC())
	  _nstep = 1;
	else
	  _nstep = 2;
	double tmax = _htime0 > 0 ? _htime0 : _htime1;
	tmax += _strawele.maxResponseTime();
	_vmax = primarywf.sampleWaveform(tmax);
	StrawDigi const& digi = digis->back();
	_tdc0 = digi.TDC(StrawDigi::zero);
	_tdc1 = digi.TDC(StrawDigi::one);
	_adc.clear();
	for(auto iadc=digi.adcWaveform().begin();iadc!=digi.adcWaveform().end();++iadc){
	  _adc.push_back(*iadc);
	}
	// mc truth information
	art::Ptr<StepPointMC> const& spmc = xpair[0]->_ihitlet->stepPointMC();
	if(!spmc.isNull()){
	  _mctime = spmc->time();
	  _mcenergy = spmc->ionizingEdep();
// compute the DOCA for this step
	  TwoLinePCA pca( straw.getMidPoint(), straw.getDirection(), 
	      spmc->position(), spmc->momentum().unit() );
	  _mcdca = pca.dca(); 
	} else {
	  _mctime = _mcenergy = _mcdca = -1000.0;
	}
	// fill the tree entry
	_sddiag->Fill();
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
DEFINE_ART_MODULE(StrawDigisFromStepPointMCs);


