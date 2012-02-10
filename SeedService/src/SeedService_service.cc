//
// Assist in the distribution of guaranteed unique seeds to all engines within a job.
//
// $Id: SeedService_service.cc,v 1.5 2012/02/10 16:26:53 gandr Exp $
// $Author: gandr $
// $Date: 2012/02/10 16:26:53 $
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
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Persistency/Provenance/ModuleDescription.h"
#include "art/Persistency/Provenance/EventID.h"

// Supporting library include files
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "cetlib/exception.h"

// C++ include files
#include <iostream>
#include <iomanip>

using namespace std;

namespace mu2e {

  namespace {
    std::string policyName[SeedService::numPolicies];

    struct NameInitializer {
      NameInitializer() {
        policyName[SeedService::unDefined]        = "unDefined";
        policyName[SeedService::autoIncrement]    = "autoIncrement";
        policyName[SeedService::preDefinedOffset] = "preDefinedOffset";
        policyName[SeedService::preDefinedSeed]   = "preDefinedSeed";
      }
    };

    NameInitializer init;
  }

  SeedService::SeedService(fhicl::ParameterSet const& pSet,
                           art::ActivityRegistry&     iRegistry  ) :

    // Initialization for all policies.
    verbosity_(pSet.get<int>("verbosity",0)),
    state_(),
    policy_(unDefined),
    pSet_(pSet),
    knownSeeds_(),
    baseSeed_(0),
    checkRange_(true),
    maxUniqueEngines_(20),

    // Initailization specific to the autoIncrement policy
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
    case preDefinedOffset: case preDefinedSeed:
      parsePreDefined();
      break;
    }

    // Register callbacks.
    iRegistry.watchPreModuleConstruction  (this, &SeedService::preModuleConstruction  );
    iRegistry.watchPostModuleConstruction (this, &SeedService::postModuleConstruction );
    iRegistry.watchPreModuleBeginRun      (this, &SeedService::preModuleBeginRun      );
    iRegistry.watchPostModuleBeginRun     (this, &SeedService::postModuleBeginRun     );
    iRegistry.watchPostEndJob             (this, &SeedService::postEndJob             );

    if ( verbosity_ > 1 ) {
      print();
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

  // Print summary information.
  void SeedService::print( ) const{

    string strCheckRange = (checkRange_) ? "true" : "false";
    mf::LogInfo log("SEEDS");
    log << "\nSummary of seeds computed by the SeedService.\n";
    log << " Policy:                       " << policyName[policy_]<< "\n";
    log << " Check range:                  " << strCheckRange     << "\n";
    log << " Maximum unique seeds per job: " << maxUniqueEngines_ << "\n";
    log << " Base Seed:                    " << baseSeed_         << "\n";
    log << " Verbosity:                    " << verbosity_        << "\n\n";

    if ( !knownSeeds_.empty() ) {
      log << " Seed Value     ModuleLabel.InstanceName\n";

      for ( map_type::const_iterator i=knownSeeds_.begin(), e=knownSeeds_.end();
            i != e; ++i ){
        log << setw(10) << i->second << "      "
            << i->first
            << "\n";
      }
    }

  }

  // Do the real work of getting and validating a seed. Not user callable.
  SeedService::seed_t SeedService::getSeed( SeedServiceHelper::EngineId const& id ){

    // Are we being called from the right place?
    ensureValidState();

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
    case autoIncrement:
      seed = currentSeed_++;
      break;
    case preDefinedOffset: case preDefinedSeed:
      seed = getPreDefined(id);
      break;
    }

    // Throw if the seed is not unique within this job.
    ensureUnique(id, seed);

    // Throw if the seed is out of range.
    ensureRange(id, seed);

    // Save the result.
    result.first->second = seed;

    return seed;
  }

  void SeedService::setPolicy(){

    string strPolicy = pSet_.get<string>("policy");
    string * iter = std::find(policyName, policyName+numPolicies, strPolicy);
    policy_ = Policy(iter - policyName);

    if ( policy_ == unDefined ){
      std::ostringstream os;
      os<< "SeedService::setPolicy(): Unrecognized policy: "
        << policyName
        << "\n Known policies are: ";

      std::copy(policyName+1, policyName+numPolicies, std::ostream_iterator<string>(os, ", "));

      throw cet::exception("SEEDS")<<os.str();
    }
  }

  void SeedService::parseAutoIncrement(){
    parseCommon();
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

    if ( checkRange_ ){
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
  void SeedService::ensureValidState( ){

    /*
     * Disable this test until the MixFilter template is upgraded to allow seeding in the c'tor.
     *
    if ( state_.state == SeedServiceHelper::ArtState::unDefined ) {
      throw cet::exception("SEEDS")
        << "SeedService: not in a module constructor or beginRun method. May not call getSeed.\n";
    }
    */
  }

  // Callbacks called by art.  Used to maintain information about state.
  void SeedService::preModuleConstruction ( art::ModuleDescription const& md){
    state_.set( SeedServiceHelper::ArtState::inConstructor, md.moduleLabel_ );
  }

  void SeedService::postModuleConstruction( art::ModuleDescription const& md){
    state_.clear();
  }

  void SeedService::preModuleBeginRun ( art::ModuleDescription const& md){
    state_.set( SeedServiceHelper::ArtState::inBeginRun, md.moduleLabel_ );
  }

  void SeedService::postModuleBeginRun( art::ModuleDescription const& md){
    state_.clear();
  }

  void SeedService::postEndJob(){
    if ( pSet_.get<bool>("endOfJobSummary",false) || verbosity_ > 1 ){
      print();
    }
  }

} // end namespace mu2e

DEFINE_ART_SERVICE(mu2e::SeedService);
