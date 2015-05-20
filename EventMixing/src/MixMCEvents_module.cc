//
// Read zero or more events from an art event-data file and mix them into the current event.
// Part of the job is to update art::Ptr objects to point into the mixed collections.
// The GenParticles, SimParticles and StatusG4 objects are always read and mixed.
// There are parameter set variables to control which of the many StepPointMCCollections
// are mixed; mixing of the PointTrajectoryCollections can also be turned on/off with a
// parameter set variable.
//
// $Id: MixMCEvents_module.cc,v 1.16 2014/03/20 18:35:10 gandr Exp $
// $Author: gandr $
// $Date: 2014/03/20 18:35:10 $
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
// 3) For additional documentation, see:
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
//      - mixStatusG4
//
// 6) The mixOp methods return a bool.  If this is true, the output product will be added to the event.
//    If it is false, the output product will not be added to the event.
//
// 7) The art::MixFilter template requires that the arguments of a mixOp methods is:
//      ( std::vector< T const*> const&, T&, art::PtrRemapper const& )
//    The second argument is not std::vector<T>&, just T&.  This might cause problems
//    if T is not a collection type.  The one case in which we encounter this is with the StatusG4
//    objects but those objects end up in the MixingSUmmary, so there is no real problem.
//
// 8) Todo:
//    When art v1_0_0 is available
//     - add extra argument to the mixOp for StatusG4.
//     - Get scale factor from the event (Done: DNB, 20 May 2015)
//

// Mu2e includes
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/MixingSummary.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/VirtualDetectorId.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Mu2eUtilities/inc/PoissonHistogramBinning.hh"
#include "SeedService/inc/SeedService.hh"

// Includes from art
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Modules/MixFilter.h"
#include "art/Framework/IO/ProductMix/MixHelper.h"
#include "art/Framework/Core/PtrRemapper.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Persistency/Common/CollectionUtilities.h"
#include "art/Utilities/InputTag.h"

// Includes from the art tool chain.
#include "cetlib/map_vector.h"

// ROOT includes
#include "TH1F.h"

// Other third party includes
#include "boost/noncopyable.hpp"
#include "CLHEP/Random/RandPoissonQ.h"

// C++ includes
#include <vector>
#include <memory>

using namespace std;

namespace mu2e {

  // Forward declare the classs we will write later in this file.e
  class MixMCEventsDetail;

  // This is the module class.
  typedef art::MixFilter<MixMCEventsDetail> MixMCEvents;

}

// Now declare the class.
class mu2e::MixMCEventsDetail : private boost::noncopyable {

public:

  MixMCEventsDetail(fhicl::ParameterSet const &p,
                    art::MixHelper &helper);

  // Called at the start of each event.
  void startEvent(const art::Event&);

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

  bool
  mixStatusG4( std::vector< mu2e::StatusG4 const *> const &in,
               mu2e::StatusG4&                            out,
               art::PtrRemapper const& );

private:

  // Run-time configurable members.
  fhicl::ParameterSet params_;

  // The number of mix-in events to choose on each event.
  double mean_;

  // Module labels of the producer that made the GenParticles and that which
  std::string genModuleLabel_;
  std::string g4ModuleLabel_;
  // label of the proton bunch intensity object
  std::string pbiLabel_;
  const ProtonBunchIntensity* pbi_;
  // Input tag for the status block.
  art::InputTag g4StatusTag_;

  // The instance names of the StepPointMCCollections to be mixed in.
  // Default is to mix all such collections.
  std::vector<mu2e::StepInstanceName> stepInstances_;

  // Enable/disable mixing of the PointTrajectoryCollections.
  bool doPointTrajectories_;
  // use proton bunch intensities or not
  bool usePBI_;

  // Only one of the following is meaningful in any instance of this class
  // If mean_ >0 then the poisson distribution is valid; else n0_ is valid.
  unique_ptr<CLHEP::RandPoissonQ> poisson_;
  int n0_;

  // Non-run-time configurable members.

  // Some counters used for diagnostics.
  int evtCount_;
  int stepCollectionCount_;

  // Histogram to record the distribution of generated events.
  TH1F* hNEvents_;

  // The number of mix-in events chosen on this event.
  size_t actual_;

  // The offsets returned from flattening the GenParticle and SimParticle collections.
  // These are needed in Ptr remapping operations.
  std::vector<size_t> genOffsets_;
  std::vector<size_t> simOffsets_;

  // New data products that will be added to the event.
  std::unique_ptr<mu2e::MixingSummary> summary_;

  // Remap all of the art::Ptr objects in one SimParticle.
  void updateSimParticle( mu2e::SimParticle& sim,
                          size_t genOffset,
                          size_t simOffset,
                          art::PtrRemapper const& remap
                          );

  // Parse the parameter set to learn which StepPointMCCollections to mix.
  std::vector<mu2e::StepInstanceName> chooseStepInstances( fhicl::ParameterSet const& pset);

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

  // Sum of the sizes of all input collections.
  // This is a candidate to be moved to a more general library.
  template <class T>
  size_t totalSize( std::vector< T const*> const& in){
    size_t sum(0);
    for ( typename std::vector<T const*>::const_iterator i=in.begin(), e=in.end();
          i !=e ; ++i ){
      sum += (*i)->size();
    }
    return sum;
  }

  // Fill output argument with size of each collection from the input argument.
  template <class T>
  void getSizes( std::vector< T const*> const& in, std::vector<size_t>& out){
    out.reserve(in.size());
    for ( typename std::vector<T const*>::const_iterator i=in.begin(), e=in.end();
          i !=e ; ++i ){
      out.push_back((*i)->size());
    }
  }

  // Variant of getSizes for use when T is a cet::map_vector; call the delta()
  // method instead of the size() method.
  template <class T>
  void getDeltas( std::vector< T const*> const& in, std::vector<size_t>& out){
    out.reserve(in.size());
    for ( typename std::vector<T const*>::const_iterator i=in.begin(), e=in.end();
          i !=e ; ++i ){
      out.push_back((*i)->delta());
    }
  }

} // end anonymous namespace

// Register mix operations and call produces<> for newly created data products
mu2e::MixMCEventsDetail::
MixMCEventsDetail(fhicl::ParameterSet const &pSet,
                  art::MixHelper &helper)
  :
  // Run-time configurable parameters
  params_(pSet.get<fhicl::ParameterSet>("detail")),
  mean_(params_.get<double>("mean")),
  genModuleLabel_(params_.get<string>("genModuleLabel")),
  g4ModuleLabel_ (params_.get<string>("g4ModuleLabel")),
  pbiLabel_ (params_.get<string>("ProtonBunchIntensityLabel")),
  g4StatusTag_   (params_.get<string>("g4StatusTag")),
  stepInstances_(chooseStepInstances(params_)),
  doPointTrajectories_(params_.get<bool>("doPointTrajectories",true)),
  usePBI_(params_.get<bool>("useProtonBunchIntensity",false)),
  poisson_(nullptr),
  n0_(0),

  // Non-run-time configurable
  evtCount_(-1),
  stepCollectionCount_(-1),
  hNEvents_(0),
  actual_(0),
  genOffsets_(),
  simOffsets_(),
  summary_(nullptr){

  // Declare new products produced directly by this class.
  helper.produces<mu2e::MixingSummary>();

  // Register MixOp operations; the callbacks are called in the order they were registered.
  helper.declareMixOp
    ( art::InputTag(genModuleLabel_,""),
      &MixMCEventsDetail::mixGenParticles, *this );

  helper.declareMixOp
    ( art::InputTag(g4ModuleLabel_,"s0"),
      &MixMCEventsDetail::mixSimParticles, *this );

  // Declare MixOps for all StepPointMCCollections.
  for ( std::vector<mu2e::StepInstanceName>::const_iterator i=stepInstances_.begin(),
          e=stepInstances_.end();  i != e; ++ i ){
    helper.declareMixOp
      ( art::InputTag(g4ModuleLabel_,i->name()),
        &MixMCEventsDetail::mixStepPointMCs, *this );
  }

  if ( doPointTrajectories_ ){
    helper.declareMixOp
      ( art::InputTag(g4ModuleLabel_,""),
        &MixMCEventsDetail::mixPointTrajectories, *this );
  }

  // When art v1_0_0 becomes available, add the extra argument so that the mixop must return false.
  if(g4StatusTag_ != art::InputTag()) {
    helper.declareMixOp
      ( g4StatusTag_,
        &MixMCEventsDetail::mixStatusG4, *this );
  }

  if ( mean_ > 0 ) {
    art::RandomNumberGenerator::base_engine_t& engine = art::ServiceHandle<art::RandomNumberGenerator>()->getEngine();
    int dummy(0);
    engine.setSeed( art::ServiceHandle<SeedService>()->getSeed(), dummy );
    poisson_ = unique_ptr<CLHEP::RandPoissonQ>( new CLHEP::RandPoissonQ(engine, mean_) );
  } else{
    n0_ = static_cast<int>(floor(std::abs(mean_)));
  }

  art::ServiceHandle<art::TFileService> tfs;
  PoissonHistogramBinning binning(mean_);
  hNEvents_ = tfs->make<TH1F>( "hNEvents", "Number of Mixed in Events", binning.nbins(), binning.xlow(), binning.xhigh() );

} // end mu2e::MixMCEventsDetail::MixMCEventsDetail

// Initialize state for each event,
void
mu2e::MixMCEventsDetail::
startEvent(const art::Event& event) {
  summary_.reset(new mu2e::MixingSummary());
  stepCollectionCount_ = -1;
  pbi_ = 0;
  // find the proton intensity for this event
  if(usePBI_){
    art::Handle<mu2e::ProtonBunchIntensity> pbiH; 
    if(event.getByLabel(pbiLabel_,pbiH))
      pbi_ = pbiH.product();
    if(pbi_ == 0)
      throw cet::exception("SIM")<<"mu2e::MixMCEvents: No ProtonBunchIntensity with label " <<  pbiLabel_ << endl;
  }
}

size_t
mu2e::MixMCEventsDetail::
nSecondaries() {

  // If the mean is positive, draw random variates from a Poisson distribution with that mean.
  // If the mean is negative, draw the same number of mix-in events on every call; that number
  // is just the absolute value of the mean, truncated to the nearest integer.
  if(mean_ > 0){
    double bmean = mean_;
    if(pbi_ != 0)
    // if we're scaling by the bunch intensity, multiply that through.  Note that the units of the
    // ProtonBunchIntensity object is # of protons hitting the target for this microbunch.
      bmean *= pbi_->intensity();
    actual_ = poisson_->fire(bmean);
  } else 
    actual_ = n0_;

  hNEvents_->Fill(actual_);

  return actual_;
}

void
mu2e::MixMCEventsDetail::
processEventIDs(art::EventIDSequence const &seq) {
  summary_->eventIDs() = seq;
}

void
mu2e::MixMCEventsDetail::
finalizeEvent(art::Event &e) {
  e.put(std::move(summary_));
}

bool
mu2e::MixMCEventsDetail::
mixGenParticles( std::vector< mu2e::GenParticleCollection const *> const& in,
                 mu2e::GenParticleCollection&                             out,
                 art::PtrRemapper const & ){

  getSizes( in, summary_->genSizes() );

  // There are no Ptr's to update; just need to flatten.
  art::flattenCollections(in, out, genOffsets_ );

  return true;
} // end mu2e::MixMCEventsDetail::mixGenParticles

bool
mu2e::MixMCEventsDetail::
mixSimParticles( std::vector< mu2e::SimParticleCollection const *> const &in,
                 mu2e::SimParticleCollection&                             out,
                 art::PtrRemapper const &remap){

  if ( in.empty()      ) return true;

  getDeltas( in, summary_->simDeltas() );

  // Flatten the input collections; does not update Ptrs.
  art::flattenCollections(in, out, simOffsets_ );

  if ( out.empty() ) return true;

  // There are art::Ptrs that need updating...

  // Tool to map sections of the output collection back to their input collection.
  Stepper<mu2e::SimParticle> inputMapper( simOffsets_, out);

  // Information about the first non-empty input container.
  StepInfo info(inputMapper.next());

  // Loop over the output collection and fix Ptrs etc.
  for ( mu2e::SimParticleCollection::iterator i=out.begin(), e=out.end();
        i !=e; ++i ){

    SimParticle& sim = i->second;

    // We have crossed the highwater mark for this input collection;
    // Get infor for the next non-empty input collection.
    if ( size_t(i->first.asInt()) >= info.high ){
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

  ++stepCollectionCount_;

  StepInstanceName instance(stepInstances_.at(stepCollectionCount_));
  getSizes( in, summary_->stepSizes(instance.id()) );

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

  getDeltas( in, summary_->pointTrajectoryDeltas() );

  size_t outSize(totalSize(in));
  out.reserve(outSize);

  typedef mu2e::PointTrajectoryCollection::key_type                         key_type;
  typedef std::vector< mu2e::PointTrajectoryCollection const *>::const_iterator Iter;

  int inputIndex(-1);
  for ( Iter i=in.begin(), e=in.end(); i !=e; ++i ){

    mu2e::PointTrajectoryCollection const& trajs(**i);
    ++inputIndex;

    for ( mu2e::PointTrajectoryCollection::const_iterator t=trajs.begin(), te=trajs.end();
          t != te; ++t ){

      key_type key(t->first);
      mu2e::PointTrajectory const& traj(t->second);

      key_type newKey = key_type(unsigned(t->first.asUint()) + simOffsets_.at(inputIndex));

      // This is redundant (unless I have made a mistake); leave it in for a while.
      if ( newKey.asInt() < out.delta() ){
        throw cet::exception("RANGE")
          << "MixMCEventsDetail::mixPointTrajectories: the key already exists: "
          << newKey
          << " Input index: " << inputIndex
          << " old key: "     << key;
      }

      // Default construct an entry in the output map_vector; hold a reference to it.
      mu2e::PointTrajectory& newtraj = out[newKey];

      // Copy the input data into the output map_vector.
      newtraj = traj;

      // Update the simId element of the PointTrajectory.
      newtraj.simId() = newKey.asInt();

    } // end loop over one input collection

  } // end loop over the collection of input collections.

  return true;
} // end mu2e::MixMCEventsDetail::mixPointTrajectories

// This does not follow the usual pattern, see note 7.
bool
mu2e::MixMCEventsDetail::
mixStatusG4( std::vector< mu2e::StatusG4 const *> const &in,
             mu2e::StatusG4&                             dummy,
             art::PtrRemapper const& ){

  // Ignore the output argument in the function signtature; instead
  // copy the input information into two data members of the MixingSummary product.
  //  - the eventStatus data member contains a copy of the input StatusG4 objects
  //  - the status data member contains a roll-up of the input StatusG4 objects
  StatusG4& status(summary_->status());
  std::vector<mu2e::StatusG4>& out(summary_->eventStatus());
  out.reserve(in.size());

  typedef std::vector< mu2e::StatusG4 const *>::const_iterator Iter;

  // There are no Ptrs to update, just copy input to output.
  for ( Iter i=in.begin(), e=in.end(); i !=e; ++i ){
    out.push_back(**i);
    status.add(**i);
  }

  // Do not add the dummy output collection to the event.
  return false;

} // end mu2e::MixMCEventsDetail::mixStatusG4

// Parse the stepInstances parameter.
std::vector<mu2e::StepInstanceName>
mu2e::MixMCEventsDetail::
chooseStepInstances( fhicl::ParameterSet const& pset){

  // Default return value.
  std::vector<mu2e::StepInstanceName>  out;

  // Build the default value: all known instance names, excluding "unknown".
  std::vector<std::string>  defaultNames;
  for ( size_t i=1; i<StepInstanceName::lastEnum; ++i ){
    defaultNames.push_back(StepInstanceName(i).name());
  }

  // Get a list of names from the parameter set; if absent, use the default.
  std::vector<std::string> names = pset.get<std::vector<std::string> >( "stepInstanceNames", defaultNames);

  // Translate strings to StepInstanceName objects.
  for ( size_t  i=0; i<names.size(); ++i ){
    StepInstanceName step(StepInstanceName::findByName( names.at(i), false));
    if ( step.id() == StepInstanceName::unknown ){
      throw cet::exception("RANGE")
        << "MixMCEventsDetail::chooseStepInstances: requested StepPointMCCollection instance name does not exist: "
        << names.at(i)
        << "\n";
    }
    out.push_back(step);
  }

  // Always do the work in a well defined order.
  sort( out.begin(), out.end() );

  // Check for duplicates.
  int nDuplicate(0);
  for ( size_t i=1; i<out.size(); ++i){
    if ( out[i] == out[i-1] ){
      ++nDuplicate;
    }
  }
  if ( nDuplicate > 0 ) {
    throw cet::exception("DUPLICATE")
      << "MixMCEventsDetail::chooseStepInstances: found duplicate StepPointMCCollection instance name(s).\n"
      << "Please fix the list stepInstanceNames parameter in your .fcl file and rerun\n";
  }

  return out;
} // end mu2e::MixMCEventsDetail::chooseStepInstances

DEFINE_ART_MODULE(mu2e::MixMCEvents);
