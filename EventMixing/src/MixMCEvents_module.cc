//
// Read zero or more events from an art event-data file and mix them into the current event.
// Part of the job is to update art::Ptr objects to point into the mixed collections.
//
// $Id: MixMCEvents_module.cc,v 1.1 2011/06/30 04:38:47 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/06/30 04:38:47 $
//
// Contact person Rob Kutschke.
//
// Notes:
//
// 1) The class name of the module class is MixMCEvents.  Therefore an fcl file
//    that uses this module must contain the parameter:
//        module_type: MixMCEvents
//
// 2) There are two parts to mixing:
//      1) The inner details of root peristency and art bookkeeping.
//      2) High level ideas about what it means to event-mix various data products.
//    The goal of the design is to separate the code into parts that require expert
//    knowledge of 1 but not 2 and other parts that require expert knowledge of 2 but not 1.
//
//    The solution is that the art developers wrote "the one true mixing module" which
//    knows all about 1 and nothing about 2.  This requires that mu2e must write a set of
//    callback functions, called by the one true mixing module, to do all of the work from part 2.
//
//    The rejected alternative was to have mu2e write the module class and provide a toolkit to
//    hide the details of 1). It was judged that even an optimally designed toolkit would still
//    require that the mu2e developer have much too much knowledge of 1).
//
//    The technology by which the callbacks are implemented is this.  The "one true mixing
//    module" is a class template.  The template argument of the class template is a so called
//    "mixing detail" class that is written by mu2e: all of the callbacks are member functions
//    of the mixing detail class.  Some of the callbacks are "just there" - the template knows to
//    look for them.  Other callbacks must be registered with the template, which is done in the
//    c'tor of the mixing detail class using calls to MixHelper::declareMixOp.
//
//    If we need to do two different kinds of mixing, we write two different mixing detail
//    classses and instantiate the template twice, once with each detail class.  Each instantiation
//    needs to be in a separate xxx_module.cc file.
//
// 3) For documentation, see:
//      test/Integration/MixAnalyzer_module.cc
//      art/Framework/Core/PtrRemapper.h
//      art/Persistency/Common/CollectionUtilities.h
//
//    The first of these is in art test suite available at:
//      https://cdcvs.fnal.gov/redmine/projects/art/repository/revisions/master/show/test/Integration
//
//    The other two are available at $ART_INC or in the code browser:
//      https://cdcvs.fnal.gov/redmine/projects/art/repository/revisions/master/show/art
//
// 4) The order in which the callbacks are called is:
//       startEvent()
//       nSecondaries();
//       processEventIDs()
//       methods registered with MixHelper::declareMixOp, in order of declaration.
//       finalizeEvent()

// Mu2e includes
#include "MCDataProducts/inc/MixingSummary.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"

// Includes from art
#include "art/Framework/Core/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Core/PtrRemapper.h"
#include "art/Persistency/Common/CollectionUtilities.h"
#include "art/Persistency/Common/PtrVector.h"
#include "art/Utilities/InputTag.h"
#include "art/Framework/Core/EngineCreator.h"

// Includes from the art tool chain.
#include "cetlib/map_vector.h"
#include "cpp0x/memory"

// Other third party includes
#include "boost/noncopyable.hpp"

#include "CLHEP/Random/RandPoissonQ.h"


using namespace std;

namespace mu2e {
  class MixMCEventsDetail;
  typedef art::MixFilter<MixMCEventsDetail> MixMCEvents;
}

class mu2e::MixMCEventsDetail : private boost::noncopyable {

public:

  MixMCEventsDetail(fhicl::ParameterSet const &p,
                    art::MixHelper &helper);

  void startEvent();

  // Return number of events to be read from the input file for one primary event.
  size_t nSecondaries();

  // Optional processEventIDs(): after the generation of the event
  // sequence, this function will be called if it exists to provide the
  // sequence of EventIDs.
  void processEventIDs(art::EventIDSequence const &seq);

  // Optional.finalizeEvent(): (eg) put bookkeping products in event. Do
  // *not* place mix products into the event: this will already have
  // been done for you.
  void finalizeEvent(art::Event &t);

  // Mixing functions. Note that they do not *have* to be member
  // functions of this detail class: they may be member functions of a
  // completely unrelated class; free functions or function objects
  // provided they (or the function object's operator()) have the
  // expected signature.

  bool aggregateGenParticleCollections( std::vector<mu2e::GenParticleCollection const *> const& in,
                                        mu2e::GenParticleCollection&                            out,
                                        art::PtrRemapper const & );

  bool
  mixSimParticles( std::vector< mu2e::SimParticleCollection const *> const &in,
                   mu2e::SimParticleCollection&                             out,
                   art::PtrRemapper const &remap);

  bool
  mixStepPointMCs( std::vector< mu2e::StepPointMCCollection const *> const &in,
                   mu2e::StepPointMCCollection&                             out,
                   art::PtrRemapper const &remap);

  bool
  mixPointTrajectories( std::vector< mu2e::PointTrajectoryCollection const *> const &in,
                        mu2e::PointTrajectoryCollection&                             out,
                        art::PtrRemapper const &remap);

private:
  double mean_;
  std::vector<size_t> doubleVectorOffsets_;
  std::vector<size_t> genOffsets_;
  std::vector<size_t> simOffsets_;
  std::auto_ptr<art::EventIDSequence> eIDs_;
  bool startEvent_called_;
  bool processEventIDs_called_;
  int  stepCollectionCount_;
  std::auto_ptr<mu2e::MixingSummary> summary_;

  // Map from order number of the 
  std::vector<size_t> whichStep_;

  CLHEP::RandPoissonQ poisson_;

  void updateSimParticle( mu2e::SimParticle& sim,
                          size_t genOffset,
                          size_t simOffset,
                          art::PtrRemapper const& remap
                          );
};

// Constructor is responsible for registering mix operations with
// MixHelper::declareMixOp() and bookkeeping products with
// MixHelperproduces().
mu2e::MixMCEventsDetail::
MixMCEventsDetail(fhicl::ParameterSet const &pSet,
                  art::MixHelper &helper)
  :
  mean_(pSet.get<double>("mean", 2)),
  doubleVectorOffsets_(),
  genOffsets_(),
  simOffsets_(),
  eIDs_(),
  startEvent_called_(false),
  processEventIDs_called_(false),
  stepCollectionCount_(0),
  summary_(0),
  whichStep_(0),
  poisson_( art::ServiceHandle<art::RandomNumberGenerator>()->getEngine(), std::abs(mean_)){

  // A mean of less than zero means to read a fixed number of events on each call.


  std::string mixProducerLabel(pSet.get<std::string>("mixProducerLabel",
                                                     "mixProducer"));
  cerr << "MixMCEventsDetail::c'tor: "
       << mean_    << " "
       << mixProducerLabel << " "
       << "\n";

  // List of StepPointMC collections in the event.
  std::vector<std::string> stepCollections;
  stepCollections.push_back("tracker");
  stepCollections.push_back("virtualdetector");
  stepCollections.push_back("stoppingtarget");
  stepCollections.push_back("CRV");
  stepCollections.push_back("calorimeter");
  stepCollections.push_back("calorimeterRO");

  // Declare new products produced directly by this class.
  // Do not declare products managed by the MixOp calls below.
  helper.produces<art::EventIDSequence>();
  helper.produces<mu2e::MixingSummary>();

  // Register MixOp operations.  Only do so for products that we
  // expect to find in these events.  MixOp operations are done
  // in the order declared here.
  helper.declareMixOp
    ( art::InputTag("generate",""),
      &MixMCEventsDetail::aggregateGenParticleCollections, *this );

  helper.declareMixOp
    ( art::InputTag("g4run",""),
      &MixMCEventsDetail::mixSimParticles, *this );

  // Declare MixOps for all StepPointMCCollections.
  for ( std::vector<std::string>::const_iterator i=stepCollections.begin(),
          e=stepCollections.end();  i != e; ++ i ){
    helper.declareMixOp
      ( art::InputTag("g4run",*i),
        &MixMCEventsDetail::mixStepPointMCs, *this );
  }

  helper.declareMixOp
    ( art::InputTag("g4run",""),
      &MixMCEventsDetail::mixPointTrajectories, *this );

  //produces<StatusG4>();

} // end mu2e::MixMCEventsDetail::MixMCEventsDetail

// Initialize state for each event,
void
mu2e::MixMCEventsDetail::
startEvent() {
  startEvent_called_ = true;
  eIDs_.reset();
  summary_.reset();
  stepCollectionCount_ = 0;
}

size_t
mu2e::MixMCEventsDetail::
nSecondaries() {

  // If the mean is negative, interpret it as a request for a fixed number of secondaries
  // on every call; that fixed number is the absolute value of the mean, rounded to the nearest
  // integer.
  static size_t nSecondaries = (mean_ < 0) ? static_cast<int>(floor(std::abs(mean_)+0.5)) : 0;

  if ( mean_ > 0 ){
    nSecondaries = poisson_.fire();
  }
  return nSecondaries;
}

void
mu2e::MixMCEventsDetail::
processEventIDs(art::EventIDSequence const &seq) {
  processEventIDs_called_ = true;
  eIDs_.reset(new art::EventIDSequence(seq));
  summary_.reset(new mu2e::MixingSummary());
  summary_->eventIDs() = seq;
  cerr << "Have my event IDs: " << eIDs_->size() << " \n";
  for ( size_t i=0; i<eIDs_->size(); ++i){
    cerr << "   "
         << i << " "
         << eIDs_->at(i) << " "
         << summary_->eventIDs().at(i) 
         << "\n";
  }
}

void
mu2e::MixMCEventsDetail::
finalizeEvent(art::Event &e) {
  cerr << "Calling finalize ... \n";
  e.put(eIDs_);
  e.put(summary_);

  assert(startEvent_called_);
  assert(processEventIDs_called_);
  startEvent_called_ = false;
  processEventIDs_called_ = false;
  cerr << "Done finalize ... " << e.event() << endl;
}

bool
mu2e::MixMCEventsDetail::
aggregateGenParticleCollections( std::vector< mu2e::GenParticleCollection const *> const& in,
                                 mu2e::GenParticleCollection&                             out,
                                 art::PtrRemapper const & ){

  cerr << "mu2e::MixMCEventsDetail::aggregateGenParticleCollections "
      << in.size()
      << "\n";
  int n(0);
  for ( std::vector< mu2e::GenParticleCollection const *>::const_iterator
          i = in.begin(), e=in.end(); i != e; ++i ){
    mu2e::GenParticleCollection const& gens = **i;
    cerr << "Next file:\n";
    for ( mu2e::GenParticleCollection::const_iterator j=gens.begin(), je=gens.end();
          j != je; ++j){
      cerr << "In:  " << n++ << " " << *j << "\n";
    }
  }
  art::flattenCollections(in, out, genOffsets_);

  cerr << "Offsets: ";
  for ( size_t i=0; i!=genOffsets_.size(); ++i){
    cerr << genOffsets_[i];
  }
  cerr << "\n";

  for ( size_t i=0; i!=out.size(); ++i){
    cerr << "out: "
        << i << " "
        << out[i]
        << "\n";
  }

  return true;
} // end mu2e::MixMCEventsDetail::aggregateGenParticleCollections

bool
mu2e::MixMCEventsDetail::
mixSimParticles( std::vector< mu2e::SimParticleCollection const *> const &in,
                 mu2e::SimParticleCollection&                             out,
                 art::PtrRemapper const &remap){
  cerr << "Mark 1 \n";

  // Nothing to do.
  if ( in.empty() ) {
    cerr << "MixMCEventDetail::mixSimParticles: Empty input container:\n";
    return true;
  }

  // Flatten the input collections; does not update Ptrs.
  art::flattenCollections(in, out, simOffsets_);

  for ( size_t i=0; i<simOffsets_.size(); ++i ){
    cerr << "Fixed Offsets: "
         << i << " "
         << in.at(i)->size() << " "
         << simOffsets_.at(i) <<  " "
         << genOffsets_.at(i)
         << "\n";
  }

  // Highwater mark that indicates one key larger than the end of the output container.
  size_t lastHighWater =  out.delta();
  cerr << "Last offset: " << lastHighWater << "\n";

  // Index into the input array and the offset arrays.
  int idx(0);

  // Highwater mark in the ouput container, marking one past the end of first input container.
  size_t highWater = (simOffsets_.size()>1 ) ? simOffsets_[1] : lastHighWater;

  // Loop over the output collection and fix Ptrs plus the internal id of each object.
  int m(0);
  for ( mu2e::SimParticleCollection::iterator i=out.begin(), e=out.end();
        i !=e; ++i ){

    SimParticle& sim = i->second;

    // We have crossed the highwater mark.  We are starting the next input container.
    // Update the highwater mark to point, within the output container, to one past the 
    // end of the next input container.
    if ( int(i->first.asInt()) >= highWater ){
      ++idx;
      highWater = ( idx<(simOffsets_.size()-1)) ? simOffsets_.at(idx+1) : lastHighWater;
    }
    cerr << "   Sim: "
        << m++ << " "
        << i->first << " "
        << sim.id() << " "
        << highWater  << " | "
        << idx  << " "
        << genOffsets_.at(idx) << " "
        << simOffsets_.at(idx) << "  |  ";
    art::Ptr<SimParticle> old = sim.parent();
    updateSimParticle( sim, genOffsets_.at(idx), simOffsets_.at(idx), remap );
    cerr << "        "
        << sim.id() << " "
        << sim.parent().id() << " "
        << sim.parent().key() << " |  "
        << old.id() << " "
        << old.key() << " : "
        << sim.genParticle().id()  << " "
        << sim.genParticle().key() << " "
        << "\n";
  }

  return true;

} // end mu2e::MixMCEventsDetail::mixSimParticles

// Update one SimParticle to deal with the flattening of the SimParticleCollections.
void
mu2e::MixMCEventsDetail::
updateSimParticle( mu2e::SimParticle& sim,
                   size_t genOffset,
                   size_t simOffset,
                   art::PtrRemapper const& remap
                   ){

  // Id of the SimParticle - must match the map_vector key.
  sim.id() = cet::map_vector_key( sim.id().asInt() + simOffset );

  // Ptr to the parent SimParticle.
  if ( sim.parent().isNonnull() ){
    sim.parent() = remap(sim.parent(), simOffset);
  }

  // Ptr to the corresponding genParticle.
  if ( sim.genParticle().isNonnull() ){
    sim.genParticle() = remap( sim.genParticle(), genOffset);
  }

  // Ptrs to all of the daughters; all will be non-null but the vector may be empty.
  std::vector<art::Ptr<SimParticle> > const& daughters= sim.daughters();
  std::vector<art::Ptr<SimParticle> > newDaughters;
  for ( size_t i=0; i != daughters.size(); ++i){
    art::Ptr<SimParticle> const& dau = remap ( daughters.at(i), simOffset );
    newDaughters.push_back( dau );
  }
  if ( !newDaughters.empty() ) sim.setDaughterPtrs( newDaughters);

} // end mu2e::MixMCEventsDetail::updateSimParticle

bool
mu2e::MixMCEventsDetail::
mixStepPointMCs( std::vector< mu2e::StepPointMCCollection const *> const &in,
                 mu2e::StepPointMCCollection&                             out,
                 art::PtrRemapper const &remap){
  cerr << "Start StepPointMCCollection: " << stepCollectionCount_++ << "\n";

  // Nothing to do.
  if ( in.empty() ) {
    cerr << "MixMCEventDetail::mixStepPointMCs: Empty input container:\n";
    return true;
  }

  // Flatten the input collections; does not update Ptrs.
  std::vector<size_t> offsets;
  art::flattenCollections(in, out, offsets);
  cerr << "Sizes: " << offsets.size() << " " << in.size() << "\n";
  cerr << "Offsets are: \n";
  for ( size_t i=0; i<offsets.size(); ++i){
    cerr <<  in.at(i)->size() << " " << offsets.at(i) << "\n";
  }

  // Highwater mark that indicates one key larger than the end of the output container.
  size_t lastHighWater =  out.size();
  cerr << "Last offset Step: " << lastHighWater << "\n";

  // Index into the input array and the offset arrays.
  int idx(0);

  // Highwater mark in the ouput container, marking one past the end of first input container.
  int highWater = (offsets.size()>1 ) ? offsets[1] : lastHighWater;
  int base(0);

  // Loop over the output collection and fix Ptrs plus the internal id of each object.
  int m(-1);
  for ( mu2e::StepPointMCCollection::iterator i=out.begin(), e=out.end();
        i !=e; ++i ){

    StepPointMC& step = *i;
    ++m;

    // We have crossed the highwater mark.  We are starting the next input container.
    // Update the highwater mark to point, within the output container, to one past the
    // end of the next input container.
    if ( m == highWater ){
      ++idx;
      base = highWater;
      highWater = ( idx<(offsets.size()-1)) ? offsets.at(idx+1) : lastHighWater;
      cerr << "Highwater: " << idx << " " << highWater << "\n";
    }
    cerr << "   Step: "
         << m << " "
         << base            << " "
         << highWater       << " "
         << idx             << " | "
         << int(step.volumeId()) << " "
         << step.time()     << " | "
         << step.simParticle().isNonnull() << " ";
    StepPointMC const& step0 = in[idx]->at(m-base);

    cerr << (step0.volumeId()) << " "
         << step0.time()       << " | ";
    if ( step.simParticle().isNonnull() ){
      art::Ptr<SimParticle> const& old = step.simParticle();
      cerr << old.id()  << " "
           << old.key() << " | " ;
      step.simParticle() = remap( step.simParticle(), simOffsets_.at(idx) );
      cerr << step.simParticle().id()  << " "
           << step.simParticle().key() << " " 
           << simOffsets_.at(idx);
    }

    cerr  << "\n";
  }

  cerr << "Done StepPointMCCollection: " << stepCollectionCount_ << "\n";
  return true;
} // end mu2e::MixMCEventsDetail::mixStepPointMCs

bool
mu2e::MixMCEventsDetail::
mixPointTrajectories( std::vector< mu2e::PointTrajectoryCollection const *> const &in,
                      mu2e::PointTrajectoryCollection&                             out,
                      art::PtrRemapper const &remap){
  cerr << "Start Point TrajectoryCollection: \n";

  std::vector<size_t> dummyOffsets;
  art::flattenCollections(in, out, dummyOffsets);
  cerr << "Sizes: " << dummyOffsets.size() << " " << in.size() << "\n";
  cerr << "DummyOffsets are: \n";
  for ( size_t i=0; i<dummyOffsets.size(); ++i){
    cerr <<  in.at(i)->size() << " "
         << dummyOffsets.at(i) << " "
         << simOffsets_.at(i)
         << "\n";
  }
  return true;
} // end mu2e::MixMCEventsDetail::mixPointTrajectories

DEFINE_ART_MODULE(mu2e::MixMCEvents);
