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
#include "Alignment/inc/AlignmentMap.hh"

namespace mu2e {
// Forward declarations
  class AlignmentSequence;

  class AlignmentService {
public:
    AlignmentService(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~AlignmentService() { delete _alignmentMap; _alignmentMap = 0;}

    // Functions registered for callbacks.
    void preBeginRun( art::Run const &run);
    void postEndJob();

    SimpleConfig const& config() const { return *_config;}

    AlignmentMap* alignmentMap() { return _alignmentMap; }
    //   AlignmentMap const * alignmentMap() const { return _alignmentMap; }

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

    // Load Alignments
    AlignmentMap* _alignmentMap;


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

