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
	fhicl::Atom<unsigned> csize{ Name("OutputCollectionSize"),
	  Comment("Estimated size of output collection"), 2000};
	fhicl::Atom<std::string> trackerSteps { Name("trackerStepPoints"),
	  Comment("Tracker StepPointMC Producer Instance name")};
      };
      using Parameters = art::EDProducer::Table<Config>;
      explicit CreateStrawGasSteps(const Parameters& conf);

    private:
      typedef std::map< std::pair<StrawId,cet::map_vector_key>, PtrStepPointMCVector > StrawSimPMap; // steps by straw, SimParticle
      void produce(art::Event& e) override;
      int _debug;
      unsigned _csize;
      bool _firstEvent;
      string _trackerStepPoints;

  };

  CreateStrawGasSteps::CreateStrawGasSteps(const Parameters& config )  : 
    art::EDProducer{config},
    _debug(config().debug()),
    _csize(config().csize()),
    _firstEvent(false),
    _trackerStepPoints(config().trackerSteps())
    {
      consumesMany<StepPointMCCollection>();
      produces <StrawGasStepCollection>();
      produces <StrawGasStepAssns>();
    }

  void CreateStrawGasSteps::produce(art::Event& event) {
  // create output
    std::unique_ptr<StrawGasStepCollection> sgsc(new StrawGasStepCollection);
    std::unique_ptr<StrawGasStepAssns> sgsa(new StrawGasStepAssns);
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
      
      StrawSimPMap stmap;

      art::Handle<StepPointMCCollection> const& handle(*ispmcc);
      StepPointMCCollection const& steps(*handle);
      nsps += steps.size();
      // Loop over the StepPointMCs in this collection and sort them
      // by straw and SimParticle
      for (size_t ispmc =0; ispmc<steps.size();++ispmc) {
	StrawId const & sid = steps[ispmc].strawId();
	// Skip dead straws, and straws that don't exist
	if (tracker.strawExists(sid)) {
	  // extract the SimParticle id
	  cet::map_vector_key tid = steps[ispmc].simParticle().get()->id();
	  // create pair
	  std::pair<StrawId,cet::map_vector_key> stpair(sid,tid);
	  // try to find it in the map
	  auto const pos = stmap.find(stpair);
	  // create ptr
	  art::Ptr<StepPointMC> spmcptr(handle,ispmc);
	  // since stmap is a Map, we loop over all steps just once and fill the map
	  // however since the map elements are pairs of paris(strawid,trackid) and a vector
	  // we need to check if a straw/track element was created already
	  if ( pos!=stmap.end() ) {
	    // the straw/track has been seen already
	    // add the ptr to the vector corresponding to this strawid and the track
	    // typedef std::vector<art::Ptr<StepPointMC> > PtrStepPointMCVector;
	    pos->second.push_back(spmcptr);
	  } else {
	    // StepPointMC seen for the first time for this straw/track
	    // create a new vector with one element;
	    // may want to reserve some based on the step limit size
	    std::vector<art::Ptr<StepPointMC>> spmcptrv(1,spmcptr);
	    stmap.emplace(stpair,spmcptrv);
	  }
	}
      }
      // now that the StepPoints are sorted by particle and straw, create the StrawGas objects and fill the collection.  Loop over the SimParticle/straw pairs
      for(auto istmap = stmap.begin(); istmap != stmap.end(); istmap++){
	auto const& spmcptrs = istmap->second;

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
	XYZVec end = Geom::toXYZVec(last->position() + last->stepLength()*last->momentum().unit());
	// average first and last time, momentum
	time = 0.5*(first->time() + last->time());
	mom = 0.5*(first->momentum().mag() + last->momentum().mag());
	// 2nd pass through steps to get width
	Hep3Vector axis;
	if(last != first)
	  axis = (last->position() - first->position()).unit();
	else
	  axis = first->momentum().unit();

	double width(0.0);
	for(auto const& spmcptr : spmcptrs)
	  width = std::max(width,(spmcptr->position()-first->position()).perp2(axis));

	// create the gas step
	StrawGasStep sgs(istmap->first.second, istmap->first.first,
		(float)eion, (float)pathlen, (float)width, (float)mom, time, 
		start, end);
	cet::map_vector_key key(istep++);
	sgsc->insert(make_pair(key,sgs));

	auto sgsp = art::Ptr<StrawGasStep>(StrawGasStepCollectionPID,sgsc->size()-1,StrawGasStepCollectionGetter);
	for(auto const& spmcptr : spmcptrs)
	  sgsa->addSingle(sgsp,spmcptr);

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
    }
    event.put(std::move(sgsc));
    event.put(std::move(sgsa));
  } // end of produce

}

DEFINE_ART_MODULE(mu2e::CreateStrawGasSteps)
