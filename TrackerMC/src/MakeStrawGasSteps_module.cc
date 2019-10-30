//
//  This module creates the StrawGasStep objects used in downstream digi simulation, using the
//  G4 StepPointMCs
//
//  Original author: David Brown (LBNL), Krzysztof Genser 19 Aug. 2019
//
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "cetlib_except/exception.h"
#include "fhiclcpp/types/Atom.h"
#include "canvas/Utilities/InputTag.h"

#include "TrackerGeom/inc/Tracker.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeometryService/inc/DetectorSystem.hh"
#include "art_root_io/TFileService.h"
#include "GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "GlobalConstantsService/inc/ParticleDataTable.hh"
#include "TrackerConditions/inc/DeadStraw.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "BTrk/BField/BField.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "MCDataProducts/inc/StrawGasStep.hh"
#include "MCDataProducts/inc/PtrStepPointMCVector.hh"
#include <utility>
// root
#include "TH1F.h"
#include "TTree.h"

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class MakeStrawGasSteps : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
	fhicl::Atom<int> debug{ Name("debugLevel"),
	  Comment("Debug Level"), 0};
	fhicl::Atom<int> diag{ Name("diagLevel"),
	  Comment("Diag Level"), 0};
	fhicl::Atom<bool> combineDeltas{ Name("CombineDeltas"),
	  Comment("Compress short delta-rays into the primary step"),true};
	fhicl::Atom<float> maxDeltaLength{ Name("MaxDeltaLength"),
	  Comment("Maximum step length for a delta ray to be compressed (mm)"),0.5};
	fhicl::Atom<float> radiusTolerance{ Name("RadiusTolerance"),
	  Comment("Tolerance to accept a point outside the straw radius (mm) "),0.5};
	fhicl::Atom<float> parabolicRotation{ Name("parabolicRotation"),
	  Comment("Maximum bending rotation to use parabolic end estimation (radians)"),0.15};
	fhicl::Atom<float> curlRotation{ Name("curlRotation"),
	  Comment("Minimum bending rotation to use curl end estimation (radians)"),2.0};
	fhicl::Atom<float> minionBG{ Name("minionBetaGamma"),
	  Comment("Minimum beta x gamma value to consider a particle minimmum-ionizing"),0.5};
	fhicl::Atom<float> curlRatio{ Name("CurlRatio"),
	  Comment("Maximum bend radius to straw radius ratio to consider a particle a curler"),1.0};
	fhicl::Atom<float> lineRatio{ Name("LineRatio"),
	  Comment("Minimum bend radius to straw radius ratio to consider a particle path a line"),10.0};
	fhicl::Atom<unsigned> csize{ Name("OutputCollectionSize"),
	  Comment("Estimated size of output collection"), 2000};
	fhicl::Atom<string> trackerSteps { Name("trackerStepPoints"),
	  Comment("Tracker StepPointMC Producer Instance name"),"tracker"};
	fhicl::Atom<string> keepDeltas { Name("KeepDeltasModule"),
	  Comment("Dont combine Deltas from this module's collection")};
	fhicl::Atom<unsigned> startSize { Name("StartSize"),
	  Comment("Starting size for straw-particle vector"),4};
	fhicl::Atom<bool> allStepsAssns{ Name("AllStepsAssns"),
	  Comment("Build the association to all the contributing StepPointMCs"),false};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit MakeStrawGasSteps(const Parameters& conf);

    private:
      void beginJob() override;
      void beginRun(art::Run& run) override;
      void produce(art::Event& e) override;
      typedef pair<StrawId,cet::map_vector_key> SSPair; // key for pair of straw, SimParticle
      typedef map< SSPair , PtrStepPointMCVector > SPSMap; // steps by straw, SimParticle
      typedef art::Handle<StepPointMCCollection> SPMCCH;
      typedef vector< SPMCCH > SPMCCHV;
      void fillMap(Tracker const& tracker,DeadStraw const& deadStraw, 
	  SPMCCH const& spmcch, SPSMap& spsmap);
      void compressDeltas(SPSMap& spsmap);
      void fillStep(PtrStepPointMCVector const& spmcptrs, Straw const& straw,
      ParticleData const* pdata, StrawGasStep& sgs, art::Ptr<StepPointMC>& spmcptr);
      void fillStepDiag(Straw const& straw, StrawGasStep const& sgs, art::Ptr<StepPointMC>const&  spmcptr);
      XYZVec endPosition(art::Ptr<StepPointMC>const& last, Straw const& straw,float charge);
      int _debug, _diag;
      bool _combineDeltas, _allAssns;
      float _maxDeltaLen, _radtol, _parrot, _curlrot;
      float _minionBG;
      float _curlmom, _linemom;
      unsigned _csize, _ssize;
      bool _firstEvent;
      string _trackerSteps, _keepDeltas;
      // cache the BField direction at the tracker center
      Hep3Vector _bdir;
      float _bnom; // BField in units of (MeV/c)/mm
      // diagnostic histograms
      TH1F *_hendrad, *_hphi;
      TTree* _sgsdiag;
      Float_t _prilen, _pridist, _elen, _erad, _epri, _esec, _partP, _brot, _width, _doca;
      vector<Float_t> _sdist;
      Int_t _npri, _nsec, _partPDG;
      ProditionsHandle<DeadStraw> _deadStraw_h;
  };

  MakeStrawGasSteps::MakeStrawGasSteps(const Parameters& config )  : 
    art::EDProducer{config},
    _debug(config().debug()),
    _diag(config().diag()),
    _combineDeltas(config().combineDeltas()),
    _allAssns(config().allStepsAssns()),
    _maxDeltaLen(config().maxDeltaLength()),
    _radtol(config().radiusTolerance()),
    _parrot(config().parabolicRotation()),
    _curlrot(config().curlRotation()),
    _minionBG(config().minionBG()),
    _curlRatio(config().curlRatio()),
    _lineRatio(config().lineRatio()),
    _csize(config().csize()),
    _ssize(config().startSize()),
    _firstEvent(false),
    _trackerSteps(config().trackerSteps()),
    _keepDeltas(config().keepDeltas())
    {
      consumesMany<StepPointMCCollection>();
      produces <StrawGasStepCollection>();
      // 2 associations: one to all the constituent StepPointMCs, the other to the first (in time) primary
      produces <StrawGasStepAssns>("Primary");
      if(_allAssns)produces <StrawGasStepAssns>("All");
    }

  void MakeStrawGasSteps::beginJob(){
    art::ServiceHandle<art::TFileService> tfs;
    if(_diag > 0){
      _hendrad = tfs->make<TH1F>("endrad","EndStep Radius",100,2.0,3.0);
      _hphi = tfs->make<TH1F>("phi","Step rotation angle",100,0.0,3.14);
      if(_diag > 1) {
	_sgsdiag=tfs->make<TTree>("sgsdiag","StrawGasStep diagnostics");
        _sgsdiag->Branch("prilen",&_prilen,"prilen/F");
        _sgsdiag->Branch("elen",&_elen,"elen/F");
        _sgsdiag->Branch("pridist",&_pridist,"pridist/F");
	_sgsdiag->Branch("brot",&_brot,"brot/F");
	_sgsdiag->Branch("erad",&_erad,"erad/F");
	_sgsdiag->Branch("epri",&_epri,"epri/F");
	_sgsdiag->Branch("esec",&_esec,"esec/F");
        _sgsdiag->Branch("width",&_width,"width/F");
        _sgsdiag->Branch("doca",&_doca,"doca/F");
        _sgsdiag->Branch("npri",&_npri,"npri/I");
        _sgsdiag->Branch("nsec",&_nsec,"nsec/I");
        _sgsdiag->Branch("partP",&_partP,"partP/F");
        _sgsdiag->Branch("partPDG",&_partPDG,"partPDG/I");
        _sgsdiag->Branch("sdist",&_sdist);
      }
    }
  }

  void MakeStrawGasSteps::beginRun( art::Run& run ){
    // get field at the center of the tracker
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    auto vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    auto bnom = bfmgr->getBField(vpoint_mu2e);
    _bdir = bnom.unit();
    // B in units of mm/MeV/c
    _bnom = bnom.mag()*BField::mmTeslaToMeVc;
    // pre-compute momentum thresholds for straight, arc, and curler
    const Tracker& tracker = *GeomHandle<Tracker>();
    const Straw& straw = tracker.getStraw(StrawId(0,0,0));
    float rstraw = straw.innerRadius();
    _curlmom = _bnom*rstraw;
    _linemom = _arcfac*_pcurl;
  }

  void MakeStrawGasSteps::produce(art::Event& event) {
    const Tracker& tracker = *GeomHandle<Tracker>();
    GlobalConstantsHandle<ParticleDataTable> pdt;
    DeadStraw const& deadStraw = _deadStraw_h.get(event.id());
    // create output
    unique_ptr<StrawGasStepCollection> sgsc(new StrawGasStepCollection);
    sgsc->reserve(_csize);
    unique_ptr<StrawGasStepAssns> sgsa_primary(new StrawGasStepAssns);
    unique_ptr<StrawGasStepAssns> sgsa_all(new StrawGasStepAssns);
    // needed for making Ptrs
    auto StrawGasStepCollectionPID = event.getProductID<StrawGasStepCollection>();
    auto StrawGasStepCollectionGetter = event.productGetter(StrawGasStepCollectionPID);
    // Get all of the tracker StepPointMC collections from the event:
    // This selector will select only data products with the given instance name.
    art::ProductInstanceNameSelector selector(_trackerSteps);
    SPMCCHV stepsHandles;
    event.getMany( selector, stepsHandles);
    //    const Tracker& tracker = *GeomHandle<Tracker>();
    // Informational message on the first event.
    if ( _firstEvent && _debug>0 ) {
      mf::LogInfo log("StrawDigiSim");
      log << "mu2e::MakeStrawGasSteps will use StepPointMCs from: \n";
      for ( SPMCCHV::const_iterator i=stepsHandles.begin(), e=stepsHandles.end();
	  i != e; ++i ){
	art::Provenance const& prov(*(i->provenance()));
	log  << "   " << prov.branchName() << "\n";
      }
      _firstEvent = false;
    }
    if(stepsHandles.empty()){
      throw cet::exception("SIM")<<"mu2e::MakeStrawGasSteps: No StepPointMC collections found for tracker" << endl;
    }
    // diagnostic counters
    unsigned nspmcs(0), nspss(0);
    // Loop over StepPointMC collections
    for( auto const& handle : stepsHandles) {
      StepPointMCCollection const& steps(*handle);
      nspmcs += steps.size();
      // see if we should compress deltas from this collection
      bool dcomp = _combineDeltas && (handle.provenance()->moduleLabel() != _keepDeltas);
      if(_debug > 1){
	if(dcomp)
	  cout << "Compressing collection " << handle.provenance()->moduleLabel() << endl;
	else
	  cout << "No compression for collection " << handle.provenance()->moduleLabel() << endl;
      }
      // Loop over the StepPointMCs in this collection and sort them by straw and SimParticle
      SPSMap spsmap; // map of step points by straw,sim particle
      fillMap(tracker,deadStraw,handle, spsmap);
      // optionally combine delta-rays that never leave the straw with their parent particle
      if(dcomp)compressDeltas(spsmap);
      nspss += spsmap.size();
      // convert the SimParticle/straw pair steps into StrawGas objects and fill the collection.  
      for(auto ispsmap = spsmap.begin(); ispsmap != spsmap.end(); ispsmap++){
	auto pid = ispsmap->first.second; // primary SimParticle
	auto const& spmcptrs = ispsmap->second;
	const Straw& straw = tracker.getStraw(ispsmap->first.first);
	ParticleData const* pdata = pdt->particle(spcptrs.front()->simParticle()->pdgId());
	StrawGasStep sgs;
	art::Ptr<StepPointMC> spmcptr;
	fillStep(spmcptrs,straw,pdata,sgs,spmcptr);
	sgsc->push_back(sgs);
	// create the Assns to the 'trigger' StepPointMC
	auto sgsp = art::Ptr<StrawGasStep>(StrawGasStepCollectionPID,sgsc->size()-1,StrawGasStepCollectionGetter);
	sgsa_primary->addSingle(sgsp,spmcptr);
	// optionall add Assns for all StepPoints, including delta-rays
	if(_allAssns){
	  for(auto const& spmcptr : spmcptrs)
	    sgsa_all->addSingle(sgsp,spmcptr);
	}
	if(_diag > 0)fillStepDiag(straw,sgs,spmcptr);
	if(_debug > 1){
	  // checks and printout
	  cout << " SGS with " << spmcptrs.size() << " steps, StrawId = " << sgs.strawId()  << " SimParticle Key = " << sgs.simParticleKey()
	    << " edep = " << eion << " pathlen = " << pathlen << " glen = " << sqrt((end-start).mag2()) << " width = " << width << " mom = " << mom
	    << " time = " << sgs.time() << endl;

	  // check if end is inside physical straw
	  const Straw& straw = tracker.getStraw(sgs.strawId());
	  static double r2 = straw.innerRadius()*straw.innerRadius();
	  Hep3Vector hend = Geom::Hep3Vec(end);
	  double rd2 = (hend-straw.getMidPoint()).perpPart(straw.getDirection()).mag2();
	  if(rd2 - r2 > 1e-5 ) cout << "End outside straw, radius " << sqrt(rd2) << endl;
	}
      } // end of pair loop
    } // end of collection loop
    if(_debug > 0){
      cout << "Total number of StrawGasSteps " << sgsc->size() << " , StepPointMCs = " << nspmcs << endl;
    }
    event.put(move(sgsc));
    event.put(move(sgsa_primary),"Primary");
    if(_allAssns) event.put(move(sgsa_all),"All");
  } // end of produce

  void MakeStrawGasSteps::fillStep(PtrStepPointMCVector const& spmcptrs, Straw const& straw,
      ParticleData const* pdata, StrawGasStep& sgs, art::Ptr<StepPointMC>& spmcptr){
    // variables we accumulate for all the StepPoints in this pair
    double eion(0.0), pathlen(0.0);
    double time(0.0), mom(0.0);
    if(_diag>1){
      _npri=_nsec=0;
      _epri=_esec=0.0;
    }
    // keep track of the first and last PRIMARY step
    art::Ptr<StepPointMC> first;
    art::Ptr<StepPointMC> last;
    // loop over all  the StepPoints for this SimParticle
    for(auto const& spmcptr : spmcptrs){
      bool primary=spmcptr->simParticle()->id() == pid;
      // update eion for all contributions
      eion += spmcptr->ionizingEdep();
      // treat primary and secondary (delta-ray) energy differently
      if(primary) {
	// primary: update path length, and entry/exit
	pathlen += spmcptr->stepLength();
	if(first.isNull() || spmcptr->time() < first->time()) first = spmcptr;
	if(last.isNull() || spmcptr->time() > last->time()) last = spmcptr;
      }
      // diagnostics
      if(_diag >1){
	if(primary){
	  _npri++;
	  _epri += spmcptr->ionizingEdep();
	} else {
	  _nsec++;
	  _esec += spmcptr->ionizingEdep();
	}
      }	    
    }
    if(first.isNull() || last.isNull())
      throw cet::exception("SIM")<<"mu2e::MakeStrawGasSteps: No first or last step" << endl;
    // for now, define the first StepPoint as the 'trigger' for this step.  Eventually
    // this might be the one closest to the wire FIXME!
    spmcptr = first;
    // Define the position at entrance and exit; note the StepPointMC position is at the start of the step, so we have to extend the last
    XYZVec start = Geom::toXYZVec(first->position());
    float charge(0.0);
    if(pdata!=0){
      charge = pdt->particle(first->simParticle()->pdgId()).ref().charge();
    }
    // determine the type of step
    StrawGasStep::StepType stype;
    setType(first,stype);
    // compute the end position and step type
    // in future we should store the end position in the StepPointMC FIXME!
    XYZVec end = endPosition(last,straw,charge,stype);
    time = first->time(); // use first time as time of this step (cluster times will be added)
    mom = 0.5*(first->momentum().mag() + last->momentum().mag());	// average first and last momentum
    // determine the width from the sigitta or curl radius
    auto pdir = first->momentum().unit();
    auto pperp = pdir.perp(_bdir);
    float bendrms = 0.5*std::min(straw.innerRadius(),mom*pperp/_bnom); // bend radius spread.  0.5 factor givs RMS of a circle
    // only sigitta perp to the wire counts
    float sint = (_bdir.cross(pdir).cross(straw.getDirection())).mag();
    static const float prms(1.0/(12.0*sqrt(5.0))); // RMS for a parabola.  This includes a factor 1/8 for the sagitta calculation too
    float sagrms = prms*sint*pathlen*pathlen*_bnom*pperp/mom;
    double width = std::min(sagrms,bendrms); // choose the smaller: different approximations work for different momenta/directions
    // create the gas step
    sgs = StrawGasStep( first->strawId(), stype,
	(float)eion,(float)pathlen, (float)width, (float)mom, time, 
	start, end);
  }

  void MakeStrawGasSteps::fillMap(Tracker const& tracker,DeadStraw const& deadStraw, 
  SPMCCH const& spmcch, SPSMap& spsmap) {
    StepPointMCCollection const& steps(*spmcch);
    for (size_t ispmc =0; ispmc<steps.size();++ispmc) {
      const auto& step = steps[ispmc];
      StrawId const & sid = step.strawId();
      Straw const& straw = tracker.getStraw(sid);
      double wpos = fabs((step.position()-straw.getMidPoint()).dot(straw.getDirection()));
      // Skip dead straws, and straws that don't exist
      //  or steps that occur in the deadened region near the end of each wire,
      // or in dead regions of the straw
      if (tracker.strawExists(sid)
	  && wpos <  straw.activeHalfLength() 
	  && deadStraw.isAlive(sid,wpos)) {
	cet::map_vector_key tid = step.simParticle().get()->id();
	// create key
	SSPair stpair(sid,tid);
	// create ptr to this step
	art::Ptr<StepPointMC> spmcptr(spmcch,ispmc);
	vector<art::Ptr<StepPointMC>> spmcptrv;
	spmcptrv.reserve(_ssize);
	spmcptrv.push_back(spmcptr);
	// check if this key exists and add it if not
	auto ssp = spsmap.emplace(stpair,spmcptrv);
	// if the key already exists, just add this Ptr to the vector
	if(!ssp.second)ssp.first->second.push_back(spmcptr);
      }
    }
  }

  void MakeStrawGasSteps::compressDeltas(SPSMap& spsmap) {
    // first, make some helper maps
    typedef map< cet::map_vector_key, StrawId > SMap; // map from key to Straw, to test for uniqueness
    typedef map< cet::map_vector_key, cet::map_vector_key> DMap; // map from delta ray to parent
    SMap smap;
    DMap dmap;
    for(auto isps = spsmap.begin(); isps != spsmap.end(); isps++ ) {
      auto sid = isps->first.first;
      auto tid = isps->first.second;
      // map key to straw
      auto sp = smap.emplace(tid,sid);
      if(!sp.second && sp.first->second != sid && sp.first->second.valid())sp.first->second = StrawId(); // Particle already seen in another straw: make invalid to avoid compressing it
    }

    // loop over particle-straw pairs looking for delta rays
    auto isps =spsmap.begin();
    while(isps != spsmap.end()){
      bool isdelta(false);
      auto& dsteps = isps->second;
      auto dkey =isps->first.second;
      // see if this particle is a delta-ray and if it's step is short
      auto pcode = dsteps.front()->simParticle()->creationCode();
      if(pcode == ProcessCode::eIoni || pcode == ProcessCode::hIoni){
	// make sure this particle doesnt have a step in any other straw
	auto ifnd = smap.find(dkey);
	if(ifnd == smap.end())
	  throw cet::exception("SIM")<<"mu2e::MakeStrawGasSteps: No SimParticle found for delta key " << dkey << endl;
	else if(ifnd->second.valid()){ // only compress delta rays without hits in any other straw
	  // add the lengths of all the steps in this straw
	  float len(0.0);
	  for(auto const& istep : dsteps)
	    len += istep->stepLength();
	  if(len < _maxDeltaLen){
	    // short delta ray. flag for combination
	    isdelta = true;
	  }
	}
      }
// if this is a delta, map it back to the primary
      if(isdelta){
	auto strawid = isps->first.first;
      // find its parent
	auto pkey = dsteps.front()->simParticle()->parentId();
	// map it so that potential daughters can map back through this particle even after compression
	dmap[dkey] = pkey;
	// find the parent. This must be recursive, as delta rays can come from delta rays (from delta rays...)
	auto jfnd = dmap.find(pkey);
	while(jfnd != dmap.end()){
	  pkey = jfnd->second;
	  jfnd = dmap.find(pkey);
	}
	// now, find the parent back in the original map
	auto ifnd = spsmap.find(make_pair(strawid,pkey));
	if(ifnd != spsmap.end()){
	  if(_debug > 1)cout << "mu2e::MakeStrawGasSteps: SimParticle found for delta parent key " << pkey << " straw " << strawid << endl;
      // move the contents to the primary
	  auto& psteps = ifnd->second;
	  psteps.insert(psteps.end(),dsteps.begin(),dsteps.end());
	  // erase the delta ray and advance the iterator
	  isps = spsmap.erase(isps);
	} else {
	// there are a very few delta rays whose parents die in the straw walls that cause StepPoints, so this is not an error.  These stay
	// as uncompressed particles.
	  if(_debug > 1)cout << "mu2e::MakeStrawGasSteps: No SimParticle found for delta parent key " << pkey << " straw " << strawid << endl;
	  isps++;
	}
      } else
      // move to the next straw/particle pair
	isps++;
    }
  }

  XYZVec MakeStrawGasSteps::endPosition(art::Ptr<StepPointMC>const& last, Straw const& straw, float charge) {
    static const double r2 = straw.innerRadius()*straw.innerRadius();
    XYZVec retval;
    // null charge has no propagation.
    if(charge == 0.0){
      retval = last->position();
    } else {
      auto momhat = last->momentum().unit();
      // test parabolic extrapolation first
      _brot = last->stepLength()*_bnom/last->momentum().mag(); // magnetic bending rotation angle
      if(_debug > 1){
	cout << "Step Length = " << last->stepLength() << " rotation angle " << _brot << endl;
      }
      auto rho = _bdir.cross(momhat);// points radially outwards for positive charge
      if(_brot < _parrot){
	// estimate end position with a parabolic trajectory
	retval  = last->position() + last->stepLength()*(momhat -(0.5*charge*_brot)*rho);
      } else {
	Hep3Vector cdir;
	if(_brot > _curlrot)
	  // curler; assume the net motion is along the BField axis.  Sign by the projection of the momentum
	  cdir = _bdir * (momhat.dot(_bdir)>0.0 ? 1.0 : -1.0);
	else
	  cdir = (momhat-(0.5*charge*_brot)*rho).unit(); // effective propagation direction
	// propagate to the straw wall
	auto pperp = (last->position()-straw.getMidPoint()).perpPart(straw.getDirection());
	double pdot = pperp.dot(cdir);
	double ppmag2 = pperp.mag2();
	double len = sqrt(pdot*pdot + r2 - ppmag2)-pdot;
	len = std::min(last->stepLength(),len);
	retval = last->position() + len*cdir;
      }
    }
    return retval; 
  }

  void MakeStrawGasSteps::fillStepDiag(Straw const& straw, StrawGasStep const& sgs,
      art::Ptr<StepPointMC>const&  spmcptr) {
    _erad = sqrt((Geom::Hep3Vec(sgs.endPosition())-straw.getMidPoint()).perpPart(straw.getDirection()).mag2());
    _hendrad->Fill(_erad);
    _hphi->Fill(_brot);
    if(_diag > 1){
      _prilen = pathlen;
      _pridist = sqrt((sgs.endPosition()-sgs.startPosition()).mag2());
      _partP = mom;
      _partPDG = spmcptr->simParticle()->pdgId();
      _elen = spmcptr->stepLength();
      _width = width;
      // compute DOCA to the wire
      TwoLinePCA poca(Geom::Hep3Vec(sgs.startPosition()),Geom::Hep3Vec(sgs.endPosition()-sgs.startPosition()),
	  straw.getMidPoint(),straw.getDirection());
      _doca = poca.dca();
      _sdist.clear();
      auto sdir = Geom::Hep3Vec(sgs.endPosition()-sgs.startPosition()).unit();
      for(auto const& spmcptr : spmcptrs){
	auto dist = ((spmcptr->position()-Geom::Hep3Vec(sgs.startPosition())).cross(sdir)).mag(); 
	_sdist.push_back(dist);
      }
      _sgsdiag->Fill();
    }
  }

  void MakeStrawGasSteps::setType(art::Ptr<StepPointMC>const& spmcptr, StrawGasStep::StepType& stype) {
  // now determine ioniztion and shape
    int itype, shape;
    if(pdata->charge() == 0.0){
      itype = StrawGasStep::StepType::neutral;
      shape = StrawGasStep::StepType::point;
    } else {
      double mom = spmcptr->momentum().mag();
      if(mom < _curlmom)
	shape = StrawGasStep::StepType::curl;
      else if(mom < _linemom)
	shape = StrawGasStep::StepType::arc;
      else
	shape = StrawGasStep::StepType::line;
      mass = pdata->mass();
      double bg = mom/mass; // beta x gamma
      if(bg > _minionBG)
	itype =StrawGasStep::StepType::minion;
      else
     	itype =StrawGasStep::StepType::highion;
    }
    stype = StrawGasStep::StepType( (StrawGasStep::StepType::shape)shape,
	(StrawGasStep::StepType::ionization)itype );
  }

}

DEFINE_ART_MODULE(mu2e::MakeStrawGasSteps)
