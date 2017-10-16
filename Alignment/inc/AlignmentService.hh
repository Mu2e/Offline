#ifndef Alignment_AlignmentService_hh
#define Alignment_AlignmentService_hh
//
// Maintain up to date alignment information and serve it to
// other services and to the modules.
//
// Original author David Norvil Brown (UofL), based
// on GeometryService by Rob Kutschke
//

// C++ include files
#include <string>
#include <memory>

// Framework include files
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art/Framework/Services/Registry/ServiceMacros.h"
#include "cetlib_except/exception.h"

#include "ConfigTools/inc/SimpleConfig.hh"
#include "boost/shared_ptr.hpp"


namespace mu2e {
// Forward declarations
  class AlignmentMap;
  class AlignmentSequence;

  class AlignmentService {
public:
    AlignmentService(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~AlignmentService();

    // Functions registered for callbacks.
    void preBeginRun( art::Run const &run);
    void postEndJob();

    SimpleConfig const& config() const { return *_config;}



private:

    // The name of the run-time configuration file.
    // Some day this will become a database key or similar.
    std::string _inputfile;

    // Control the behaviour of messages from the SimpleConfig object holding
    // the geometry parameters.
    bool _allowReplacement;
    bool _messageOnReplacement;
    bool _messageOnDefault;
    int  _configStatsVerbosity;

    // Print final config file after all replacements.
    bool _printConfig;

    // The object that parses run-time configuration file.
    std::unique_ptr<SimpleConfig> _config;

    // Load G4 geometry options
    std::unique_ptr<AlignmentMap> _alignmentMap;

    // Check the configuration.
    //    void checkConfig();

    // template <typename DET> friend class GeomHandle;

    // template <class DET>
    // DET* getElement()
    // {
    //   if(_run_count==0)
    //     throw cet::exception("GEOM")
    //       << "Cannot get detectors from an unconfigured geometry service.\n"
    //       << "You've attempted to a get an element before the first run\n";

    //   // to use this generic way requires a map of names (typeid?) to
    //   // abstract elements.
    //   // find the detector element requested
    //   std::string name = typeid(DET).name();
    //   DetMap::iterator it(_detectors.find(name));
    //   if(it==_detectors.end())
    //     throw cet::exception("GEOM")
    //       << "Failed to retrieve detector element of type " << name << "\n";

    //   // this must succeed or there is something terribly wrong
    //   DET* d = dynamic_cast<DET*>(it->second.get());

    //   if(d==0)
    //     throw cet::exception("GEOM")
    //       << "Failed to convert found detector " << name
    //       << " to its correct type.  There is a serious problem.\n";

    //   return d;
    // }


    // Keep a count of how many runs we have seen.
    int _run_count;

    // This is not copyable or assignable - private and unimplemented.
    AlignmentService const& operator=(AlignmentService const& rhs);
    AlignmentService(AlignmentService const& rhs);

    //    // Don't need to expose definition of private template in header
    //    template <typename DET> void addDetector(std::unique_ptr<DET> d);
    //    template <typename DETALIAS, typename DET> void addDetectorAliasToBaseClass(std::unique_ptr<DET> d);

  };

}

DECLARE_ART_SERVICE(mu2e::AlignmentService, LEGACY)
#endif // Alignment_AlignmentService_hh

