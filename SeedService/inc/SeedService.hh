#ifndef SeedService_SeedService_hh
#define SeedService_SeedService_hh
//
// An art service to assist in the distribution of guaranteed unique seeds to all engines within an art job.
//
// $Id: SeedService.hh,v 1.1 2012/02/04 00:12:48 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/02/04 00:12:48 $
//
// Contact person Rob Kutschke
//
// The following instructions presume that the user is familiar with the information on:
//    http://mu2e.fnal.gov/atwork/computing/Random.shtml#background
//    http://mu2e.fnal.gov/atwork/computing/Random.shtml#user
//
// This class is configured from a fhcil parameter set:
//
//    SeedService : {
//       policy         : "autoIncrement"  // Required: Other legal value is "preDefined"
//       baseSeed       : 0                // Required: An integer >= 0.
//       checkRange     : true             // Optional: legal values true or false; defaults to true
//       maxSeedsPerJob : 20               // Required iff checkRange is true.
//
//       verbosity       : 0               // Optional: default=0, no informational printout
//       endOfJobSummary : false           // Optional: print list of all managed seeds at end of job.
//    }
//
// The policy parameter tells the service to choose one of two alogrithms, "autoIncrement"
// of "preDefined".  We anticipate that other algorithms may be added in the future.
// If the value of the policy parameter is not one of the known policies, the code will
// throw an exception.
//
// If the policy is defined as "autoIncrement" the above fhcil fragment shows all
// run time configurable items. Additional parameters for the "preDefined" policy
// will be discussed later.
//
// Module code requests a seed by making one of the following two calls:
//
//    art::ServiceHandle<SeedService>->getSeed();
//    art::ServiceHandle<SeedService>->getSeed("instanceName");
//
// It is the callers responsibility to use the appropriate form. The call to getSeed must be
// done in either the constructor or the beginRun method of a module.  When a call is made
// to getSeed, the SeedService is aware of which module instance is making the call; each
// module instance is identified by its module label.
//
// When getSeed is called with a particular module label and instance name, it computes
// a seed value, saves it and returns it.  If there is a subsequent call to getSeed with
// the same module label and instance name, the service will return the saved value of
// the seed.  The following text will use the phrase "unique calls to getSeed"; two
// calls with the same module label and instance names are not considered unique.
//
// If the policy is set to autoIncrement, the seed is set to baseSeed+offset, where
// on the first unique call to getSeed the offset is set to 0; on the second unique
// call to getSeed it is set to 1, and so on.
//
// If the policy is set to preDedefined then, when getSeed is called, the service will look
// into the parameter set to find a defined offset for the specified module label
// and instance name.  The returned value of the seed will be baseSeed+offset.
//
// The fhicl grammar to specify the offsets takes two forms.  If no instance name
// is given, the offset is given by:
//
//   moduleLabel : offset
//
// When a module has multiple instances, the offsets are given by:
//
//   moduleLabel : {
//      instanceName1 : offset1
//      instanceName2 : offset2
//   }
//
// The SeedService does several additional checks.
//
// If one (module label, instance name) has the same seed as another (module label, instance name),
// the service will throw an exception.
//
// If the checkRange parameter is set to true, and if an offset is generated with a value outside
// the range ( 0<= offset < maxSeedsPerJob-1 ) then the service will throw an exception.
//
// It is the responsibility of the user to ensure that the baseSeed and maxSeedsPerJob are chosen
// it a way that ensures the required level of uniqueness of seeds.  The example grid jobs have
// a single point of maintenance to achieve this: the user must specify the starting run number
// for each grid submission.
//

// Some helper classes.
#include "SeedService/inc/ArtState.hh"
#include "SeedService/inc/EngineId.hh"

// From art and its tool chain.
#include "art/Framework/Services/Optional/RandomNumberGenerator.h"
#include "fhiclcpp/ParameterSet.h"

#include <string>
#include <map>

// Forward declarations
namespace art {
  class ActivityRegistry;
  class ModuleDescription;
  class Run;
  class SubRun;
}

namespace mu2e {

  class SeedService {
  public:

    typedef art::RandomNumberGenerator::seed_t seed_t;

    enum Policy { unDefined, autoIncrement, preDefined};

    SeedService(const fhicl::ParameterSet&, art::ActivityRegistry&);

    // Accept compiler written d'tor.  Not copyable or assignable.

    // Return the seed value for this module label (instance name).
    seed_t getSeed();
    seed_t getSeed( std::string const& instanceName );

    // Print known (EngineId,seed) pairs.
    void print() const;

    // Call backs that will be called by art.
    void preModuleConstruction (art::ModuleDescription const& md);
    void postModuleConstruction(art::ModuleDescription const& md);
    void preModuleBeginRun     (art::ModuleDescription const& md);
    void postModuleBeginRun    (art::ModuleDescription const& md);
    void postEndJob();

  private:

    // This class is not copyable or assignable: these methods are not implemented.
    SeedService(SeedService const& rhs);
    SeedService const& operator=(SeedService const& rhs);

    // Control the level of information messages.
    int verbosity_;

    // Helper to track state of art. It is legal to call getSeed only in certain states.
    SeedServiceHelper::ArtState state_;

    // Which of the supported policies to use?
    Policy policy_;

    // Copy of the parameter set passed in at c'tor time.
    fhicl::ParameterSet pSet_;

    // List of seeds already computed.  Seeds are unique per EngineID (not per call to getSeed).
    typedef std::map<SeedServiceHelper::EngineId,seed_t> map_type;
    map_type knownSeeds_;

    // Information used by all policies.
    seed_t baseSeed_;
    bool   checkRange_;
    seed_t maxSeedsPerJob_;

    // Information used by the autoIncrement policy
    seed_t currentSeed_;

    // Main logic for computing and validating a seed.
    seed_t getSeed( SeedServiceHelper::EngineId const& );

    // Helper functions for all policies
    void setPolicy       ( );
    void ensureValidState( );
    void ensureRange     ( SeedServiceHelper::EngineId const& id, seed_t seed );
    void ensureUnique    ( SeedServiceHelper::EngineId const& id, seed_t seed );
    void parseCommon     ();

    // Helper functions for the autoIncrement policy.
    void parseAutoIncrement();

    // Helper functions for the preDefined policy.
    void   parsePreDefined   ();
    seed_t getPreDefinedSeed ( SeedServiceHelper::EngineId const& id );

  };

}

#endif /* SeedService_SeedService_hh */
