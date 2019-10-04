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
#include "BFieldGeom/inc/BFieldManager.hh"
#include "BTrk/BField/BField.hh"

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/MCRelationship.hh"
#include "MCDataProducts/inc/StrawGasStep.hh"
#include "MCDataProducts/inc/PtrStepPointMCVector.hh"
#include <utility>

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
	fhicl::Atom<bool> combineDeltas{ Name("CombineDeltas"),
	  Comment("Compress short delta-rays into the primary step"),true};
	fhicl::Atom<float> maxDeltaLength{ Name("MaxDeltaLength"),
	  Comment("Maximum step length for a delta ray to be compressed"),0.5}; //mm
	fhicl::Atom<float> radiusTolerance{ Name("RadiusTolerance"),
	  Comment("Tolerance to accept a point outside the straw radius"),0.5}; //mm
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
      void beginRun(art::Run& run) override;
      void produce(art::Event& e) override;
      typedef pair<StrawId,cet::map_vector_key> SSPair; // key for pair of straw, SimParticle
      typedef map< SSPair , PtrStepPointMCVector > SPSMap; // steps by straw, SimParticle
      void compressDeltas(SPSMap& spsmap);
      XYZVec endPosition(art::Ptr<StepPointMC>const& last, Straw const& straw);
      int _debug;
      bool _combineDeltas, _allAssns;
      float _maxDeltaLen, _radtol;
      unsigned _csize, _ssize;
      bool _firstEvent;
      string _trackerSteps, _keepDeltas;
      // cache the BField direction at the tracker center
      Hep3Vector _bdir, _bnom;
  };

  MakeStrawGasSteps::MakeStrawGasSteps(const Parameters& config )  : 
    art::EDProducer{config},
    _debug(config().debug()),
    _combineDeltas(config().combineDeltas()),
    _allAssns(config().allStepsAssns()),
    _maxDeltaLen(config().maxDeltaLength()),
    _radtol(config().radiusTolerance()),
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

  void MakeStrawGasSteps::beginRun( art::Run& run ){
    // get field at the center of the tracker
    GeomHandle<BFieldManager> bfmgr;
    GeomHandle<DetectorSystem> det;
    Hep3Vector vpoint_mu2e = det->toMu2e(Hep3Vector(0.0,0.0,0.0));
    _bnom = bfmgr->getBField(vpoint_mu2e);
    if ( _bnom.mag() < 1.0e-4 ){
      _bdir = Hep3Vector(0.0,0.0,1.0);
    } else {
      _bdir = _bnom.unit();
    } 
  }

  void MakeStrawGasSteps::produce(art::Event& event) {
    const Tracker& tracker = *GeomHandle<Tracker>();
    // create output
    unique_ptr<StrawGasStepCollection> sgsc(new StrawGasStepCollection);
    sgsc->reserve(_csize);
    unique_ptr<StrawGasStepAssns> sgsa_primary(new StrawGasStepAssns);
    unique_ptr<StrawGasStepAssns> sgsa_all(new StrawGasStepAssns);
    // needed for making Ptrs
    auto StrawGasStepCollectionPID = event.getProductID<StrawGasStepCollection>();
    auto StrawGasStepCollectionGetter = event.productGetter(StrawGasStepCollectionPID);
    // Get all of the tracker StepPointMC collections from the event:
    typedef vector< art::Handle<StepPointMCCollection> > HandleVector;
    // This selector will select only data products with the given instance name.
    art::ProductInstanceNameSelector selector(_trackerSteps);
    HandleVector stepsHandles;
    event.getMany( selector, stepsHandles);
    //    const Tracker& tracker = *GeomHandle<Tracker>();
    // Informational message on the first event.
    if ( _firstEvent && _debug>0 ) {
      mf::LogInfo log("StrawDigiSim");
      log << "mu2e::MakeStrawGasSteps will use StepPointMCs from: \n";
      for ( HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end();
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
      SPSMap spsmap; // map of step points by straw,sim particle
      StepPointMCCollection const& steps(*handle);
      nspmcs += steps.size();
      // see if we should compress deltas from this collection
      bool dcomp = _combineDeltas && handle.provenance()->moduleLabel() != _keepDeltas;
      if(_debug > 1 && !dcomp)cout << "No compression for collection " << handle.provenance()->moduleLabel() << endl;
      // Loop over the StepPointMCs in this collection and sort them by straw and SimParticle
      for (size_t ispmc =0; ispmc<steps.size();++ispmc) {
	const auto& step = steps[ispmc];
	StrawId const & sid = step.strawId();
	// Skip dead straws, and straws that don't exist
	if (tracker.strawExists(sid)) {
	  cet::map_vector_key tid = step.simParticle().get()->id();
	  // create key
	  pair<StrawId,cet::map_vector_key> stpair(sid,tid);
	  // create ptr to this step
	  art::Ptr<StepPointMC> spmcptr(handle,ispmc);
	  vector<art::Ptr<StepPointMC>> spmcptrv;
	  spmcptrv.reserve(_ssize);
	  spmcptrv.push_back(spmcptr);
	  // check if this key exists and add it if not
	  auto ssp = spsmap.emplace(stpair,spmcptrv);
	  // if the key already exists, just add this Ptr to the vector
	  if(!ssp.second)ssp.first->second.push_back(spmcptr);
	}
      }
      // optionally combine delta-rays that never leave the straw with their parent particle
      if(dcomp)compressDeltas(spsmap);
      nspss += spsmap.size();
      // convert the SimParticle/straw pair steps into StrawGas objects and fill the collection.  
      for(auto ispsmap = spsmap.begin(); ispsmap != spsmap.end(); ispsmap++){
	auto const& spmcptrs = ispsmap->second;
	const Straw& straw = tracker.getStraw(ispsmap->first.first);
	// variables we accumulate for all the StepPoints in this pair
	double eion(0.0), pathlen(0.0);
	double time(0.0), mom(0.0);
	// keep track of the first and last PRIMARY step
	art::Ptr<StepPointMC> first(spmcptrs.front());
	art::Ptr<StepPointMC> last(spmcptrs.front());
	// loop over all  the StepPoints for this SimParticle
	auto pid = first->trackId();
	for(auto const& spmcptr : spmcptrs){
	// update eion for all contributions
	  eion += spmcptr->ionizingEdep();
	// treat primary and secondary (delta-ray) energy differently
	  if(spmcptr->simParticle()->id() == pid) {
	  // primary: update path length, and entry/exit
	    pathlen += spmcptr->stepLength();
	    if(spmcptr->time() < first->time()) first = spmcptr;
	    if(spmcptr->time() > last->time()) last = spmcptr;
	  }
	}
	// Define the position at entrance and exit; note the StepPointMC position is at the start of the step, so we have to extend the last
	XYZVec start = Geom::toXYZVec(first->position());
	XYZVec end = endPosition(last,straw);
	time = first->time(); // use first time as time of this step (cluster times will be added)
	mom = 0.5*(first->momentum().mag() + last->momentum().mag());	// average first and last momentum
	// determine the width from the triangle calculation
	double width(0.0);
	if(spmcptrs.size() > 1)width = sqrt( pathlen*pathlen - (end-start).mag2())/2.0;

	// create the gas step
	StrawGasStep sgs(ispsmap->first.second, ispsmap->first.first,
		(float)eion,(float)pathlen, (float)width, (float)mom, time, 
		start, end);
	sgsc->push_back(sgs);
	// create the Assns
	auto sgsp = art::Ptr<StrawGasStep>(StrawGasStepCollectionPID,sgsc->size()-1,StrawGasStepCollectionGetter);
	sgsa_primary->addSingle(sgsp,first);
	// add Assns for all StepPoints, including delta-rays, for diagnostics
	if(_allAssns){
	  for(auto const& spmcptr : spmcptrs)
	    sgsa_all->addSingle(sgsp,spmcptr);
	}

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

  XYZVec MakeStrawGasSteps::endPosition(art::Ptr<StepPointMC>const& last, Straw const& straw) {
// test linear extrapolation first
    auto momhat = last->momentum().unit();
    auto end = last->position() + last->stepLength()*momhat;
    // check position against straw boundary
    static const double r2 = straw.innerRadius()*straw.innerRadius();
    double rd2 = (end-straw.getMidPoint()).perpPart(straw.getDirection()).mag2();
    if(rd2 > r2) {
    // outside the straw: adjust the end position
      Hep3Vector mperp;
      if(rd2 - r2 < _radtol){
	// close enogh; interpolate linearly along the momentum direction;
	mperp = momhat.perpPart(straw.getDirection());
      } else {
	// curler; assume the net motion is along the BField axis.  Sign by the projection of the momentum
	mperp = _bdir;
	if(momhat.dot(_bdir)<0.0) mperp *= -1.0;
	momhat = mperp;
      }
      auto pperp = (last->position()-straw.getMidPoint()).perpPart(straw.getDirection());
      double pmdot = pperp.dot(mperp);
      double ppmag2 = pperp.mag2();
      double mpmag2 = mperp.mag2();
      double len = (sqrt(pmdot*pmdot + mpmag2*(r2 - ppmag2))-pmdot)/mpmag2;
      end = last->position() + len*momhat;
      if(_debug > 1){
	double newr2 =(end-straw.getMidPoint()).perpPart(straw.getDirection()).mag2();
	cout <<"Forced in straw, old length = " << last->stepLength() << " new length " << len << " new r2 " << newr2 << endl;
      }

    } 
    return Geom::toXYZVec(end);
  }
}

DEFINE_ART_MODULE(mu2e::MakeStrawGasSteps)
