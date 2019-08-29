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

#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawGasStep.hh"
#include "MCDataProducts/inc/PtrStepPointMCVector.hh"
#include <utility>

using namespace std;
using CLHEP::Hep3Vector;
namespace mu2e {

  class CreateStrawGasSteps : public art::EDProducer {
    public:
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      struct Config {
	fhicl::Atom<int> debug{ Name("debugLevel"),
	  Comment("Debug Level"), 0};
	fhicl::Atom<bool> compDeltas{ Name("CompressDeltas"),
	  Comment("Compress short delta-rays into the primary step")};
	fhicl::Atom<float> maxDeltaLength{ Name("MaxDeltaLength"),
	  Comment("Maximum step length for a delta ray to be compressed"),0.5};
	fhicl::Atom<unsigned> csize{ Name("OutputCollectionSize"),
	  Comment("Estimated size of output collection"), 2000};
	fhicl::Atom<string> trackerSteps { Name("trackerStepPoints"),
	  Comment("Tracker StepPointMC Producer Instance name")};
	  fhicl::Atom<unsigned> startSize { Name("StartSize"),
	  Comment("Starting size for straw-particle vector"),1};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit CreateStrawGasSteps(const Parameters& conf);

    private:
      typedef map< pair<StrawId,cet::map_vector_key>, PtrStepPointMCVector > StrawSimPMap; // steps by straw, SimParticle
      typedef map< cet::map_vector_key, StrawId> SimPStrawMap; // map from SimParticle to Straw
      void produce(art::Event& e) override;
      void extractDeltas(SimPStrawMap const& spmap,StrawSimPMap& sspmap, StrawSimPMap& deltas);
      void addDeltas(StrawGasStep& step,StrawSimPMap& deltas);
      void endPosition(art::Ptr<StepPointMC>const& last, Straw const& straw, XYZVec& lastpos);
      int _debug;
      bool _compressDeltas;
      float _maxDeltaLen;
      unsigned _csize, _ssize;
      bool _firstEvent;
      string _trackerStepPoints;

  };

  CreateStrawGasSteps::CreateStrawGasSteps(const Parameters& config )  : 
    art::EDProducer{config},
    _debug(config().debug()),
    _compressDeltas(config().compDeltas()),
    _maxDeltaLen(config().maxDeltaLength()),
    _csize(config().csize()),
    _ssize(config().startSize()),
    _firstEvent(false),
    _trackerStepPoints(config().trackerSteps())
    {
      consumesMany<StepPointMCCollection>();
      produces <StrawGasStepCollection>();
      produces <StrawGasStepAssns>();
    }

  void CreateStrawGasSteps::produce(art::Event& event) {
  // create output
    unique_ptr<StrawGasStepCollection> sgsc(new StrawGasStepCollection);
    unique_ptr<StrawGasStepAssns> sgsa(new StrawGasStepAssns);
    // needed for making Ptrs
    auto StrawGasStepCollectionPID = event.getProductID<StrawGasStepCollection>();
    auto StrawGasStepCollectionGetter = event.productGetter(StrawGasStepCollectionPID);
    // extend the collection
    sgsc->reserve(_csize);
    // Get all of the tracker StepPointMC collections from the event:
    typedef vector< art::Handle<StepPointMCCollection> > HandleVector;
    // This selector will select only data products with the given instance name.
    art::ProductInstanceNameSelector selector(_trackerStepPoints);
    HandleVector stepsHandles;
    event.getMany( selector, stepsHandles);
    const Tracker& tracker = *GeomHandle<Tracker>();
    // Informational message on the first event.
    if ( _firstEvent && _debug>0 ) {
      mf::LogInfo log("StrawDigiSim");
      log << "mu2e::CreateStrawGasSteps will use StepPointMCs from: \n";
      for ( HandleVector::const_iterator i=stepsHandles.begin(), e=stepsHandles.end();
	  i != e; ++i ){
	art::Provenance const& prov(*(i->provenance()));
	log  << "   " << prov.branchName() << "\n";
      }
      _firstEvent = false;
    }
    if(stepsHandles.empty()){
      throw cet::exception("SIM")<<"mu2e::CreateStrawGasSteps: No StepPointMC collections found for tracker" << endl;
    }
    unsigned istep(0);
    unsigned nsps(0);
    // Loop over StepPointMC collections
    for ( HandleVector::const_iterator ispmcc=stepsHandles.begin(), espmcc=stepsHandles.end();ispmcc != espmcc; ++ispmcc ){
      StrawSimPMap sspmap, deltas;
      SimPStrawMap spmap;
      art::Handle<StepPointMCCollection> const& handle(*ispmcc);
      StepPointMCCollection const& steps(*handle);
      nsps += steps.size();
      // Loop over the StepPointMCs in this collection and sort them by straw and SimParticle
      for (size_t ispmc =0; ispmc<steps.size();++ispmc) {
	StrawId const & sid = steps[ispmc].strawId();
	// Skip dead straws, and straws that don't exist
	if (tracker.strawExists(sid)) {
	  // extract the SimParticle id
	  cet::map_vector_key tid = steps[ispmc].simParticle().get()->id();
	  // check if this particle has already been seen in a different straw
	  auto sp = spmap.emplace(tid,sid);
	  if(!sp.second && sp.first->second != sid) 
	    if(sp.first->second.valid())sp.first->second = StrawId(); // Already seen in another straw: make invalid to avoid compressing it
	  // create key
	  pair<StrawId,cet::map_vector_key> stpair(sid,tid);
	  // create ptr to this step
	  art::Ptr<StepPointMC> spmcptr(handle,ispmc);
	  vector<art::Ptr<StepPointMC>> spmcptrv(_ssize,spmcptr);
	  // check if this key exists and add it if not
	  auto ssp = sspmap.emplace(stpair,spmcptrv);
	  // if the key already exists, just add this Ptr to the vector
	  if(!ssp.second)ssp.first->second.push_back(spmcptr);
	}
      }
      // optionally compress out delta-rays that never leave the straw
      
      if(_compressDeltas)extractDeltas(spmap,sspmap,deltas);
      // create the StrawGas objects and fill the collection.  Loop over the SimParticle/straw pairs
      for(auto isspmap = sspmap.begin(); isspmap != sspmap.end(); isspmap++){
	auto const& spmcptrs = isspmap->second;
	const Straw& straw = tracker.getStraw(isspmap->first.first);

	// variables we accumulate for all the StepPoints in this pair
	double eion(0.0), pathlen(0.0);
	double time(0.0), mom(0.0);
	// keep track of the first and last step
	art::Ptr<StepPointMC> first(spmcptrs.front());
	art::Ptr<StepPointMC> last(spmcptrs.front());
	// loop over all  the StepPoints for this SimParticle
	for(auto const& spmcptr : spmcptrs){
	  eion += spmcptr->ionizingEdep();
	  pathlen += spmcptr->stepLength();
	  if(spmcptr->time() < first->time()) first = spmcptr;
	  if(spmcptr->time() > last->time()) last = spmcptr;
	}
	// Define the position at entrance and exit; note the StepPointMC position is at the start of the step, so we have to extend the last
	XYZVec start = Geom::toXYZVec(first->position());
	XYZVec end;
	endPosition(last,straw,end);
	time = first->time(); // use first time as time of this step (cluster times will be added)
	mom = 0.5*(first->momentum().mag() + last->momentum().mag());	// average first and last momentum
	// 2nd pass through steps to get width
	Hep3Vector axis;
	if(last != first)
	  axis = (last->position() - first->position()).unit();
	else
	  axis = first->momentum().unit();

	double width(0.0);
	for(auto const& spmcptr : spmcptrs)
	  width = max(width,(spmcptr->position()-first->position()).perp2(axis));

	// create the gas step
	StrawGasStep sgs(isspmap->first.second, isspmap->first.first,
		(float)eion, (float)pathlen, (float)width, (float)mom, time, 
		start, end);
	cet::map_vector_key key(istep++);
	// Add in the energy from the delta-rays
	if(_compressDeltas) addDeltas(sgs,deltas);
	sgsc->insert(make_pair(key,sgs));

	// create the Assns
	auto sgsp = art::Ptr<StrawGasStep>(StrawGasStepCollectionPID,sgsc->size()-1,StrawGasStepCollectionGetter);
	for(auto const& spmcptr : spmcptrs)
	  sgsa->addSingle(sgsp,spmcptr);
	  // add delta-rays FIXME!

	if(_debug > 1){
	  // checks and printout
	  cout << " SGS with " << spmcptrs.size() << " steps, StrawId = " << sgs.strawId()  << " SimParticle Key = " << sgs.simParticleKey()
	  << " edep = " << eion << " pathlen = " << pathlen << " width = " << width << " mom = " << mom
	  << " time = " << sgs.time() << endl;
//	  << " start = " << start << endl
//	  << " end   = " << end << endl;

	  // check if end is inside physical straw
	  const Straw& straw = tracker.getStraw(sgs.strawId());
	  double r2 = straw.innerRadius()*straw.innerRadius();
	  Hep3Vector hend = Geom::Hep3Vec(end);
          double rd2 = (hend-straw.getMidPoint()).perpPart(straw.getDirection()).mag2();
	  if(rd2 - r2 > 1e-5 ) cout << "End outside straw " << rd2 << endl;
 
	
	}
      } // end of pair loop
    } // end of collection loop
    if(_debug > 0){
      cout << "Total number of StrawGasSteps " << sgsc->size() << " , StepPointMCs = " << nsps << endl;
//      << " SimParticles before " << nsimold << endl;
    }
    event.put(move(sgsc));
    event.put(move(sgsa));
  } // end of produce

  void CreateStrawGasSteps::extractDeltas(SimPStrawMap const& spmap,StrawSimPMap& sspmap, StrawSimPMap& deltas) {
    auto issp = sspmap.begin();
    while(issp != sspmap.end()){
      bool isdelta(false);
    // see if this particle is a delta-ray and if it's step is short
      auto pcode = issp->second.front()->simParticle()->creationCode();
      if(pcode == ProcessCode::eIoni || pcode == ProcessCode::hIoni){
	float len(0.0);
	for(auto const& istep : issp->second)
	  len += istep->stepLength();
	if(len < _maxDeltaLen){
// make sure this particle didn't create hits in any other straw
	  isdelta = spmap.find(issp->first.second)->second == issp->first.first;
	}
      }
// if this is a delta, move the object to the delta-ray list, and erase it from the primary map
      if(isdelta){
	deltas.emplace(*issp);
	sspmap.erase(issp++);
      } else
	++issp;
    }
  }

  void CreateStrawGasSteps::addDeltas(StrawGasStep& step,StrawSimPMap& deltas) {
    //FIXME!

  }

  void CreateStrawGasSteps::endPosition(art::Ptr<StepPointMC>const& last, Straw const& straw, XYZVec& lastpos) {
    lastpos = Geom::toXYZVec(last->position() + last->stepLength()*last->momentum().unit());
    // check straw boundary
//    double r2 = straw.innerRadius()*straw.innerRadius();
//    Hep3Vector hend = Geom::Hep3Vec(end);
//    double rd2 = (hend-straw.getMidPoint()).perpPart(straw.getDirection()).mag2();
    // FIXME!
  }
}

  DEFINE_ART_MODULE(mu2e::CreateStrawGasSteps)
