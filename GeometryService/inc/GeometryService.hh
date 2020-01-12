#ifndef GeometryService_GeometryService_hh
#define GeometryService_GeometryService_hh
//
// Maintain up to date geometry information and serve it to
// other services and to the modules.
//
// Original author Rob Kutschke
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
#include "Mu2eInterfaces/inc/Detector.hh"
#include "boost/shared_ptr.hpp"

// FIXME: Make a backdoor to geom svc to instantiate detector by hand. - call from G4_Module::beginRun.
// right after geom initialize.
//

namespace mu2e {

// Forward declarations
  class Target;
  class G4;
  class Mu2eG4Study;
  class Mu2eHall;
  class G4GeometryOptions;

  class GeometryService {
public:
    GeometryService(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~GeometryService();

    // Functions registered for callbacks.
    void preBeginRun( art::Run const &run);
    void postEndJob();

    SimpleConfig const& config() const { return *_config;}

    fhicl::ParameterSet const& pset() const { return _pset; }

    G4GeometryOptions const * geomOptions() const { return _g4GeomOptions.get(); }
    G4GeometryOptions       * geomOptions()       { return _g4GeomOptions.get(); }

    template <class DET>
    bool hasElement()
    {
      if(_run_count==0)
        throw cet::exception("GEOM")
          << "Cannot get detectors from an unconfigured geometry service.\n"
          << "You've attempted to a get an element before the first run\n";

      // to use this generic way requires a map of names (typeid?) to
      // abstract elements.
      // find the detector element requested
      std::string name = typeid(DET).name();
      DetMap::iterator it(_detectors.find(name));

      return !(it==_detectors.end());

    }

    bool isStandardMu2eDetector() const { return standardMu2eDetector_; }

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

    const fhicl::ParameterSet          _pset;

    // Load G4 geometry options
    std::unique_ptr<G4GeometryOptions> _g4GeomOptions;

    // Check the configuration.
    void checkConfig();
    void checkTrackerConfig();

    typedef boost::shared_ptr<Detector> DetectorPtr;
    typedef std::map<std::string,DetectorPtr> DetMap;

    template <typename DET> friend class GeomHandle;

    template <class DET>
    DET* getElement()
    {
      if(_run_count==0)
        throw cet::exception("GEOM")
          << "Cannot get detectors from an unconfigured geometry service.\n"
          << "You've attempted to a get an element before the first run\n";

      // to use this generic way requires a map of names (typeid?) to
      // abstract elements.
      // find the detector element requested
      std::string name = typeid(DET).name();
      DetMap::iterator it(_detectors.find(name));
      if(it==_detectors.end())
        throw cet::exception("GEOM")
          << "Failed to retrieve detector element of type " << name << "\n";

      // this must succeed or there is something terribly wrong
      DET* d = dynamic_cast<DET*>(it->second.get());

      if(d==0)
        throw cet::exception("GEOM")
          << "Failed to convert found detector " << name
          << " to its correct type.  There is a serious problem.\n";

      return d;
    }

    // is this the standard Mu2e detector?
    bool standardMu2eDetector_;

    // All of the detectors that we know about.
    DetMap _detectors;

    // Keep a count of how many runs we have seen.
    int _run_count;

    // This is not copyable or assignable - private and unimplemented.
    GeometryService const& operator=(GeometryService const& rhs);
    GeometryService(GeometryService const& rhs);

    // Don't need to expose definition of private template in header
    template <typename DET> void addDetector(std::unique_ptr<DET> d);
    template <typename DETALIAS, typename DET> void addDetectorAliasToBaseClass(std::unique_ptr<DET> d);

    // Some information that is provided through the GeometryService
    // should only be used inside GEANT jobs.  The following method is
    // used by G4 to make this info available.
    friend class G4;
    friend class Mu2eG4;
    friend class Mu2eG4MT; 
    friend class Mu2eG4Study;
    void addWorldG4(const Mu2eHall&);

  };

}

DECLARE_ART_SERVICE(mu2e::GeometryService, SHARED)
#endif /* GeometryService_GeometryService_hh */
