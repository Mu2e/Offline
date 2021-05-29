//
// Assist in the distribution of guaranteed unique seeds to all engines within a job.
//
//
// Contact person Rob Kutschke
//

// Todo:
// 1) Move end of job printout to a print(ostream& ) method.
// 2) Switch to turn printout on and off
// 3) Add verbose switch.


// Mu2e include files
#include "SeedService/inc/SeedService.hh"

// Art include files
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"
#include "art/Persistency/Provenance/ModuleContext.h"
#include "art/Persistency/Provenance/ModuleDescription.h"
#include "canvas/Persistency/Provenance/EventID.h"

// Supporting library include files
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib_except/exception.h"

// C++ include files
#include <iostream>
#include <iomanip>

using namespace std;

namespace mu2e {

  const std::vector<std::string>& SeedService::policyNames() {
    static std::vector<std::string> names;
    if(names.empty()) {
      const char *cnames[] = {
#define X(x) #x,
        SEED_SERVICE_POLICIES
#undef X
      };
      names = std::vector<std::string>(cnames, cnames + sizeof(cnames)/sizeof(cnames[0]));
    }

    return names;
  }

  SeedService::SeedService(fhicl::ParameterSet const& pSet,
                           art::ActivityRegistry&     iRegistry  ) :

    // Initialization for all policies.
    verbosity_(pSet.get<int>("verbosity",0)),
    state_(),
    policy_(unDefined),
    pSet_(pSet),
    knownSeeds_(),
    baseSeed_(8),
    checkRange_(true),
    maxUniqueEngines_(20),

    // Initailization specific to the autoIncrement and linearMapping policies
    currentSeed_(0){

    // Throw if policy is not recognized.
    setPolicy();

    // Finish parsing the parameter set, as required by the selected policy.
    switch(policy_) {
    default:
      throw cet::exception("SEEDS")<< "SeedService(): Internal error: unknown policy_ value\n";
    case autoIncrement:
      parseAutoIncrement();
      break;
    case linearMapping:
      parseLinearMapping();
      break;
    case preDefinedOffset: case preDefinedSeed:
      parsePreDefined();
      break;
    }

    // Register callbacks.
    iRegistry.sPreModuleConstruction.watch  (this, &SeedService::preModuleConstruction  );
    iRegistry.sPostModuleConstruction.watch (this, &SeedService::postModuleConstruction );
    iRegistry.sPreModuleBeginRun.watch      (this, &SeedService::preModuleBeginRun      );
    iRegistry.sPostModuleBeginRun.watch     (this, &SeedService::postModuleBeginRun     );
    iRegistry.sPostEndJob.watch             (this, &SeedService::postEndJob             );

    if ( verbosity_ > 0 ) {
      print(std::cout);
    }

  }

  // User callable entry, without instance name.
  SeedService::seed_t SeedService::getSeed(){
    SeedServiceHelper::EngineId id(state_.currentModuleLabel);
    return getSeed(id);
  }

  // User callable entry, with instance name.
  SeedService::seed_t SeedService::getSeed( std::string const& instanceName ){
    SeedServiceHelper::EngineId id(state_.currentModuleLabel,instanceName);
    return getSeed(id);
  }

  // User callable entry for InputSource
  SeedService::seed_t SeedService::getInputSourceSeed(){
    SeedServiceHelper::EngineId id("_source_");
    return getSeed(id);
  }

  // Print summary information.
  void SeedService::print( ) const{
    mf::LogInfo log("SEEDS");
    print(log);
  }

  // Do the real work of getting and validating a seed. Not user callable.
  SeedService::seed_t SeedService::getSeed( SeedServiceHelper::EngineId const& id ){

    // Are we being called from the right place?
    ensureValidState(id);

    // Check for an already computed seed.
    std::pair<map_type::iterator,bool> result = knownSeeds_.insert(pair<SeedServiceHelper::EngineId,seed_t>(id,0));
    if ( !result.second ){
      return result.first->second;
    }

    // Compute the seed.
    long seed(0);
    switch(policy_) {
    default:
      throw cet::exception("SEEDS")<< "getSeed(): Internal SeedService: error: unknown policy_ value\n";
    case autoIncrement: case linearMapping:
      seed = currentSeed_++;
      break;
    case preDefinedOffset: case preDefinedSeed:
      seed = getPreDefined(id);
      break;
    }

    if(policy_ != preDefinedSeed) {
      // Throw if the seed is not unique within this job.
      ensureUnique(id, seed);

      // Throw if the seed is out of range.
      ensureRange(id, seed);
    }

    // Save the result.
    result.first->second = seed;

    return seed;
  }

  void SeedService::setPolicy(){

    const string strPolicy = pSet_.get<string>("policy");
    std::vector<std::string>::const_iterator iter = std::find(policyNames().begin(), policyNames().end(), strPolicy);
    if(iter != policyNames().end()) {
      policy_ = Policy(std::distance(policyNames().begin(), iter));
    }

    if ( policy_ == unDefined ){
      std::ostringstream os;
      os<< "SeedService::setPolicy(): Unrecognized policy: "
        << strPolicy
        << "\n Known policies are: ";

      std::copy(policyNames().begin(), policyNames().end(), std::ostream_iterator<string>(os, ", "));

      throw cet::exception("SEEDS")<<os.str();
    }
  }

  void SeedService::parseAutoIncrement(){
    parseCommon();
    currentSeed_  = baseSeed_;
  }

  void SeedService::parseLinearMapping(){
    parseCommon();
    baseSeed_ *= maxUniqueEngines_;
    currentSeed_  = baseSeed_;
  }

  void SeedService::parsePreDefined(){
    parseCommon();
  }


  // Parse elements that are common to all policies
  void SeedService::parseCommon(){

    if ( !pSet_.get_if_present("baseSeed",baseSeed_) ){
      throw cet::exception("SEEDS")
        << "SeedService: was unable to find the parameter baseSeed in the input parameter set.\n";
    }
    checkRange_ = pSet_.get<bool>("checkRange",true);

    if ( checkRange_ || (policy_ == linearMapping) ){
      maxUniqueEngines_ = pSet_.get<seed_t>("maxUniqueEngines");
    }

  }

  // This method handles both the preDefinedOffset and preDefinedSeed cases
  SeedService::seed_t SeedService::getPreDefined( SeedServiceHelper::EngineId const& id ){

    // Case 1: no instance name.
    if ( !id.instanceDefined ){
      seed_t offset(0);
      if ( !pSet_.get_if_present( id.moduleLabel, offset) ){
        throw cet::exception("SEEDS")
          << "SeedService: is in preDefinedOffset mode and was unable to find the seed offset for"
          << id.moduleLabel << "\n";
      }
      return (policy_ == preDefinedOffset) ? baseSeed_ + offset : offset;
    }

    // Case 2: instance name is given.
    fhicl::ParameterSet subSet;
    if ( !pSet_.get_if_present( id.moduleLabel, subSet) ){
      throw cet::exception("SEEDS")
        << "SeedService: is in preDefinedOffset mode and was unable to find the seed block for module: "
        << id.moduleLabel
        << " and instance name: "
        << id.instanceName
        << "\n";
    }

    seed_t offset(0);
    if ( !subSet.get_if_present( id.instanceName, offset) ){
      throw cet::exception("SEEDS")
        << "SeedService: is in preDefinedOffset mode and was unable to find the seed offset for module: "
        << id.moduleLabel
        << " and instance name: "
        << id.instanceName
        << "\n";
    }

    return (policy_ == preDefinedOffset) ? baseSeed_ + offset : offset;
  }

  // Throw if the seed has already been used.
  void SeedService::ensureUnique( SeedServiceHelper::EngineId const& id,
                                  SeedService::seed_t seed){

    for ( map_type::const_iterator i=knownSeeds_.begin(), e=knownSeeds_.end();
          i != e; ++i ){

      // Do not compare to self.
      if ( i->first == id ) continue;

      if ( i->second == seed ){
        throw cet::exception("SEEDS")
          << "SeedService::ensureUnique offset: " << seed-baseSeed_
          <<", seed: "<<seed
          << ". Already used by module.instance: " << i->first << "\n"
          << "May not be reused by module.instance: " << id << "\n";
      }
    }

  }

  // Throw if the seed exceeds the maximum allowed seed value for this job.
  void SeedService::ensureRange( SeedServiceHelper::EngineId const& id,
                                 SeedService::seed_t                seed){

    if ( !checkRange_ ) return;

    seed_t offset = seed - baseSeed_;
    if ( seed >= baseSeed_ + maxUniqueEngines_ ){
      throw cet::exception("SEEDS")
        << "SeedService:ensureRange for engine: "
        << id << " the seed offset is: "
        << offset
        << ".\nAllowed seed offsets are in the range 0....(N-1) where N is: "
        << maxUniqueEngines_ <<"\n";
    }
  }

  // getSeed may only be called from a c'tor or from a beginRun method. In all other cases, throw.
  void SeedService::ensureValidState(SeedServiceHelper::EngineId const& id){

    if( (state_.state == SeedServiceHelper::ArtState::unDefined) &&
        // A kludge to allow use from a Source constructor.
        // Note that the underscore can never occur in a real module label.
        (id.moduleLabel != "_source_")
        ) {
      throw cet::exception("SEEDS")
        << "SeedService: not in a module constructor or beginRun method. May not call getSeed.\n";
    }

  }

  // Callbacks called by art.  Used to maintain information about state.
  void SeedService::preModuleConstruction ( art::ModuleDescription const& md){
    state_.set( SeedServiceHelper::ArtState::inConstructor, md.moduleLabel() );
  }

  void SeedService::postModuleConstruction( art::ModuleDescription const&){
    state_.clear();
  }

  void SeedService::preModuleBeginRun ( art::ModuleContext const& mc){
    state_.set( SeedServiceHelper::ArtState::inBeginRun, mc.moduleLabel() );
  }

  void SeedService::postModuleBeginRun( art::ModuleContext const&){
    state_.clear();
  }

  void SeedService::postEndJob(){
    if(verbosity_ > 0) {
      print(std::cout);
    }
    else if(pSet_.get<bool>("endOfJobSummary",false)) {
      print(); // framework logger decides whether and where it shows up
    }
  }

} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::SeedService);
