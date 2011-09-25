//
// Read zero or more events from an art event-data file and mix them into the current event.
// Part of the job is to update art::Ptr objects to point into the mixed collections.
//
// $Id: MixMCEvents_module.cc,v 1.2 2011/09/25 18:39:09 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/09/25 18:39:09 $
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
//    If we want to mix DIOs from one file and EjectedProtons from another file, then the .fcl
//    file needs to instantiate two instances of the same mixing module, giving them different
//    input files and, presumably, different Poisson means.
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
//       - startEvent()
//       - nSecondaries();
//       - processEventIDs()
//       - methods registered with MixHelper::declareMixOp, in order of declaration.
//       - finalizeEvent()
//
// 5) In this module, the order of declareMixOp methods is:
//      - mixGenParticles
//      - mixSimParticles
//      - mixStepPointMCs, possibly called many times
//      - mixPointTrajectories
//      - coming soon: accessing all of the G4status objects. Not 100% sure how this will be done?
//
// 6) The mixOp methods return a bool.  If this is true, the output product will be added to the event.
//    If it is false, the output product will not be added to the event.
//

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
#include "art/Utilities/InputTag.h"
//#include "art/Framework/Core/EngineCreator.h"

// Includes from the art tool chain.
#include "cetlib/map_vector.h"
#include "cpp0x/memory"

// Other third party includes
#include "boost/noncopyable.hpp"

#include "CLHEP/Random/RandPoissonQ.h"


#include <vector>

using namespace std;

namespace mu2e {
  class MixMCEventsDetail;

  // This is the module class.
  typedef art::MixFilter<MixMCEventsDetail> MixMCEvents;

}

class mu2e::MixMCEventsDetail : private boost::noncopyable {

public:

  MixMCEventsDetail(fhicl::ParameterSet const &p,
                    art::MixHelper &helper);

  // Called at the start of each event.
  void startEvent();

  // Return number of events to be read from the input file for this primary event.
  size_t nSecondaries();

  // Save the event Ids taken from the mix-in input file.
  void processEventIDs(art::EventIDSequence const &seq);

  // The last method called for each event.
  void finalizeEvent(art::Event &t);

  // Mixing functions:
  bool mixGenParticles( std::vector<mu2e::GenParticleCollection const *> const& in,
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

  // Run-time configurable members.

  // The number of mix-in events to choose on each event.
  double mean_;

  // The instance names of the StepPointMCCollections to be mixed in.
  std::vector<std::string> stepInstanceNames_;

  // Non-run-time configurable members.

  // Some counters used for diagnostics.
  int evtCount_;
  int stepCollectionCount_;

  // The number of mix-in events chosen on this event.
  size_t actual_;

  // The offsets returned from flattening the GenParticle and SimParticle collections.
  std::vector<size_t> genOffsets_;
  std::vector<size_t> simOffsets_;

  // New data products that will be added to the event.
  std::auto_ptr<art::EventIDSequence> eIDs_;
  std::auto_ptr<mu2e::MixingSummary> summary_;

  // Remap all of the art::Ptr objects in one SimParticle.
  void updateSimParticle( mu2e::SimParticle& sim,
                          size_t genOffset,
                          size_t simOffset,
                          art::PtrRemapper const& remap
                          );

  // A hack that will be needed until the mixing template is fixed.
  // Skip events if the actual number of mix-in events is different than the
  // requested number.
  bool hack( size_t );

};

// Some helper classes used only within this file.
namespace {

  // When one or more input collections is mixed into one output
  // collection, this class describes the range in the output collection that
  // corresponds to one input collection.
  struct StepInfo {
    StepInfo( size_t aindex, size_t asize, size_t alow, size_t ahigh):
      index(aindex),
      size(asize),
      low(alow),
      high(ahigh){
    }

    size_t index;    // Index identifying which element of the input collection this is.
    size_t size;     // "Size" of this input element; this has different meanings
                     //   - if the collection is an std::vector this is the size of the input collection.
                     //   - if the collection is a cet::map_vector this is the difference between the
                     //     first and last keys in the input collection.

                     // Range of indices/keys in the output collection the correspond
                     // to the elements of this input collection:
    size_t low;      //   low end of the range
    size_t high;     //   one more than the high end of the range.
  };

  inline std::ostream& operator<<(std::ostream& ost,
                                  const StepInfo& id ){
    ost << "( "
        << id.index     << " : "
        << id.size      << " : "
        << id.low      << " : "
        << id.high 
        << " )";
    return ost;
  }

  // When looping through the output collection, this class is used to help
  // keep track of where we are among the input collections.
  template<class T>
  class Stepper {
  public:

    // Two different c'tors for different collection types.
    // They different in the call to .size() or .delta().
    Stepper( std::vector<size_t > const& offsets,
             std::vector<T> const&       out
               ):
      next_(0),
      info_(){

      completeConstruction ( offsets, out.size() );
    }

    Stepper( std::vector<size_t > const& offsets,
             cet::map_vector<T> const&   out
               ):
      next_(0),
      info_(){

      completeConstruction ( offsets, out.delta() );
    }

    // Return the information about the next input collection with non-zero contents.
    StepInfo const& next() {

      for ( size_t i=next_; i!=info_.size(); ++i){
        StepInfo const& stepInfo = info_[i];
        if ( stepInfo.size > 0 ){
          next_ = i+1;
          return info_.at(i);
        }
      }

      // The user of this class should ensure that we never drop through the loop.
      throw cet::exception("RANGE")
        << "Stepper::next : error navigating mixed collection to find Ptrs.";
    }


    StepInfo const& info( size_t i) const{
        return info_.at(i);
    }

  private:

    size_t next_;
    std::vector<StepInfo> info_;

    void completeConstruction( std::vector<size_t > const& offsets, size_t lastHighWater ){
      info_.reserve(offsets.size());

      for ( size_t i=0; i!=offsets.size()-1; ++i){
        int s = offsets.at(i+1)-offsets.at(i);
        info_.push_back( StepInfo(i,s,offsets.at(i),offsets.at(i+1)));
      }
      int s = lastHighWater-offsets.back();
      info_.push_back( StepInfo(offsets.size()-1, s, offsets.back(), lastHighWater ));
    }

  }; // end class Stepper

} // end anonymous namespace

// Register mix operations and call produces<> for newly created data products
mu2e::MixMCEventsDetail::
MixMCEventsDetail(fhicl::ParameterSet const &pSet,
                  art::MixHelper &helper)
  :

  // Run-time configurable parameters
  mean_(pSet.get<double>("mean", 2.)),
  stepInstanceNames_(),                   //  Still to be made run-time-configurable.

  // Non-run-time configurable
  evtCount_(-1),
  stepCollectionCount_(0),
  actual_(0),
  genOffsets_(),
  simOffsets_(),
  eIDs_(),
  summary_(0){

  // List of StepPointMC collections in the event.
  stepInstanceNames_.push_back("tracker");
  stepInstanceNames_.push_back("virtualdetector");
  stepInstanceNames_.push_back("stoppingtarget");
  stepInstanceNames_.push_back("CRV");
  stepInstanceNames_.push_back("calorimeter");
  stepInstanceNames_.push_back("calorimeterRO");

  // Declare new products produced directly by this class.
  helper.produces<art::EventIDSequence>();
  helper.produces<mu2e::MixingSummary>();

  // Register MixOp operations; the callbacks are called in the order they were registered.
  helper.declareMixOp
    ( art::InputTag("generate",""),
      &MixMCEventsDetail::mixGenParticles, *this );

  helper.declareMixOp
    ( art::InputTag("g4run",""),
      &MixMCEventsDetail::mixSimParticles, *this );

  // Declare MixOps for all StepPointMCCollections.
  for ( std::vector<std::string>::const_iterator i=stepInstanceNames_.begin(),
          e=stepInstanceNames_.end();  i != e; ++ i ){
    helper.declareMixOp
      ( art::InputTag("g4run",*i),
        &MixMCEventsDetail::mixStepPointMCs, *this );
  }

  helper.declareMixOp
    ( art::InputTag("g4run",""),
      &MixMCEventsDetail::mixPointTrajectories, *this );

  //produces<StatusG4>();

} // end mu2e::MixMCEventsDetail::MixMCEventsDetail

// If we should skip this event because nSecondaries returned 0, then
// return true.  This will go away when we fix the bug in the mixing
// class template.
bool mu2e::MixMCEventsDetail::hack ( size_t n ){
  if ( n != actual_ ){
    cerr << "Skipping mixOp at event: "
         << evtCount_  << " : " 
         << n          << " "
         << actual_
         << endl;
    return true;
  }
  return false;
}

// Initialize state for each event,
void
mu2e::MixMCEventsDetail::
startEvent() {
  eIDs_.reset();
  summary_.reset();
  stepCollectionCount_ = 0;
}

size_t
mu2e::MixMCEventsDetail::
nSecondaries() {

  // If the mean is positive, draw random variates from a Poisson distribution with that mean.
  // If the mean is negative, draw the same number of mix-in events on every call; that number
  // is just the absolute value of the mean, truncated to the nearest integer.

  int n(0);
  if ( mean_ > 0 ){
    static CLHEP::RandPoissonQ poisson( art::ServiceHandle<art::RandomNumberGenerator>()->getEngine(), std::abs(mean_));
    n = poisson.fire();
  }else{
    static size_t n0 = (mean_ < 0) ? static_cast<int>(floor(std::abs(mean_))) : 0;
    n = n0;
  }

  actual_ = n;

  return n;
}

void
mu2e::MixMCEventsDetail::
processEventIDs(art::EventIDSequence const &seq) {
  eIDs_.reset(new art::EventIDSequence(seq));
  summary_.reset(new mu2e::MixingSummary());
  summary_->eventIDs() = seq;
}

void
mu2e::MixMCEventsDetail::
finalizeEvent(art::Event &e) {
  e.put(eIDs_);
  e.put(summary_);
}

bool
mu2e::MixMCEventsDetail::
mixGenParticles( std::vector< mu2e::GenParticleCollection const *> const& in,
                 mu2e::GenParticleCollection&                             out,
                 art::PtrRemapper const & ){


  if ( hack(in.size()) ) return true;

  // There are no Ptr's to update; just need to flatten.
  art::flattenCollections(in, out, genOffsets_);

  return true;
} // end mu2e::MixMCEventsDetail::mixGenParticles

bool
mu2e::MixMCEventsDetail::
mixSimParticles( std::vector< mu2e::SimParticleCollection const *> const &in,
                 mu2e::SimParticleCollection&                             out,
                 art::PtrRemapper const &remap){

  if ( hack(in.size()) ) return true;
  if ( in.empty()      ) return true;

  // Flatten the input collections; does not update Ptrs.
  art::flattenCollections(in, out, simOffsets_);

  if ( out.empty() ) return true;

  // There are art::Ptrs that need updating...

  // Tool to map sections of the output collection back to their input collection.
  Stepper<mu2e::SimParticle> inputMapper(simOffsets_, out);

  // Information about the first non-empty input container.
  StepInfo info(inputMapper.next());

  // Loop over the output collection and fix Ptrs etc.
  for ( mu2e::SimParticleCollection::iterator i=out.begin(), e=out.end();
        i !=e; ++i ){

    SimParticle& sim = i->second;

    // We have crossed the highwater mark for this input collection; 
    // Get infor for the next non-empty input collection.
    if ( int(i->first.asInt()) >= info.high ){
      info = inputMapper.next();
    }

    // Update Ptrs etc.
    updateSimParticle( sim, genOffsets_.at(info.index), simOffsets_.at(info.index), remap );
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
  // Redo this to avoid the copy.  Needs an update to SimParticle.hh.
  std::vector<art::Ptr<SimParticle> > const& daughters = sim.daughters();
  if ( !daughters.empty() ) {
    std::vector<art::Ptr<SimParticle> > newDaughters;
    newDaughters.reserve(daughters.size());
    for ( size_t i=0; i != daughters.size(); ++i){
      art::Ptr<SimParticle> const& dau = remap ( daughters.at(i), simOffset );
      newDaughters.push_back( dau );
    }
    sim.setDaughterPtrs( newDaughters);
  }

} // end mu2e::MixMCEventsDetail::updateSimParticle

bool
mu2e::MixMCEventsDetail::
mixStepPointMCs( std::vector< mu2e::StepPointMCCollection const *> const &in,
                 mu2e::StepPointMCCollection&                             out,
                 art::PtrRemapper const &remap){

  if ( hack(in.size()) ) return true;
  if ( in.empty()      ) return true;

  // Flatten the input collections; does not update Ptrs.
  std::vector<size_t> offsets;
  art::flattenCollections(in, out, offsets);

  if ( out.empty() ) return true;

  // There are art::Ptrs that need updating ...

  // Tool to map sections of the output collection back to their input collection.
  Stepper<mu2e::StepPointMC> nn(offsets, out);

  // Information about the first non-empty input container.
  StepInfo info(nn.next());

  // Loop over the output collection and fix Ptrs etc.
  int m(-1);
  for ( mu2e::StepPointMCCollection::iterator i=out.begin(), e=out.end();
        i !=e; ++i ){

    StepPointMC& step = *i;

    // We have crossed the highwater mark for this input collection; 
    // Get info for the next non-empty input collection.
    if ( ++m == int(info.high) ){
      info = nn.next();
    }

    if ( step.simParticle().isNonnull() ){
      step.simParticle() = remap( step.simParticle(), simOffsets_.at(info.index) );
    }

  }

  return true;
} // end mu2e::MixMCEventsDetail::mixStepPointMCs

bool
mu2e::MixMCEventsDetail::
mixPointTrajectories( std::vector< mu2e::PointTrajectoryCollection const *> const &in,
                      mu2e::PointTrajectoryCollection&                             out,
                      art::PtrRemapper const &remap){

  if  ( hack(in.size()) ) return true;

  // Build the flattened collection in a temporary std::map.
  // This will be inserted into the output collection in a separate step.
  // We will need to extend the interface of cet::map_vector in order to build the
  // the output in place and avoid the extra copies.

  typedef mu2e::PointTrajectoryCollection::key_type                             key_type;
  typedef std::vector< mu2e::PointTrajectoryCollection const *>::const_iterator Iter;
  typedef std::map<key_type,mu2e::PointTrajectory>                              tmp_type;

  tmp_type tmp;
  int inputIndex(-1);
  for ( Iter i=in.begin(), e=in.end(); i !=e; ++i ){

    mu2e::PointTrajectoryCollection const& trajs(**i);
    ++inputIndex;
    
    for ( mu2e::PointTrajectoryCollection::const_iterator t=trajs.begin(), te=trajs.end();
          t != te; ++t ){

      key_type key(t->first);
      mu2e::PointTrajectory const& traj(t->second);
    
      key_type newKey = key_type(unsigned(t->first.asUint()) + simOffsets_.at(inputIndex));

      std::pair<tmp_type::iterator,bool> ret=
        tmp.insert(std::pair<key_type,mu2e::PointTrajectory>(newKey,traj));

      // The key previously existed, there is a mistake somewhere.
      if ( ! ret.second ){
        throw cet::exception("RANGE")
          << "MixMCEventsDetail::mixPointTrajectories: the key was already in the temporary map: " 
          << newKey
          << " Input index: " << inputIndex 
          << " old key: " << key;
      }

      // Update the simId data member in the transient map.  
      // This data member will become an art::Ptr<SimParticle> some day.
      mu2e::PointTrajectory& newtraj = ret.first->second;
      newtraj.simId() = newKey.asInt();

    } // end loop over one input collection

  } // end loop over the collection of input collections.

  // Copy the temporary map into the output.
  out.insert( tmp.begin(), tmp.end());

  return true;
} // end mu2e::MixMCEventsDetail::mixPointTrajectories

DEFINE_ART_MODULE(mu2e::MixMCEvents);
