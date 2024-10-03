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
#include "art/Framework/Services/Registry/ServiceTable.h"
#include "art/Framework/Services/Registry/ServiceDeclarationMacros.h"
#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Table.h"
#include "cetlib_except/exception.h"

#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"
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
    // Validated fhcil configuration
    struct SimulatedDetector {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> tool_type{Name("tool_type"),"Mu2e"};
    };

    struct Config {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;
      fhicl::Atom<std::string> inputFile{Name("inputFile")};
      fhicl::Atom<std::string> bFieldFile{Name("bFieldFile")};
      fhicl::Atom<bool>   allowReplacement{Name("allowReplacement"),true};
      fhicl::Atom<bool>   messageOnReplacement{Name("messageOnReplacement"),false};
      fhicl::Atom<bool>   messageOnDefault{Name("messageOnDefault"),false};
      fhicl::Atom<int>    configStatsVerbosity{Name("configStatsVerbosity"),false};
      fhicl::Atom<bool>   printConfig{Name("printConfig"),false};
      fhicl::Atom<bool>   printConfigTopLevel{Name("printConfigTopLevel"),false};
      fhicl::Table<SimulatedDetector> simulatedDetector{Name("simulatedDetector")};
    };

    using Parameters= art::ServiceTable<Config>;

    GeometryService(const Parameters&, art::ActivityRegistry&);
    ~GeometryService();

    // This is not copyable or assignable.
    GeometryService(GeometryService const& ) = delete;
    GeometryService(GeometryService&&      ) = delete;
    GeometryService& operator=(GeometryService const& ) = delete;
    GeometryService& operator=(GeometryService&&      ) = delete;

    // Functions registered for callbacks.
    void preBeginRun( art::Run const &run);
    void postEndJob();

    SimpleConfig const& config() const { return *_config;}

    fhicl::ParameterSet const& simulatedDetector() const { return _simulatedDetector; }

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

    bool isStandardMu2eDetector() const { return _standardMu2eDetector; }

private:

    // The name of the run-time configuration file.
    // Some day this will become a database key or similar.
    std::string _inputfile;

    // The name of the file that contains the configuration of the B fields.
    std::string _bFieldFile;

    // Control the behaviour of messages from the SimpleConfig object holding
    // the geometry parameters. These affect both SimpleConfig objects.
    bool _allowReplacement;
    bool _messageOnReplacement;
    bool _messageOnDefault;
    int  _configStatsVerbosity;

    // Print final config file after all replacements.  These affect both SimpleConfig objects.
    bool _printConfig;
    bool _printTopLevel;

    // The objects that parse the run-time configuration files.
    std::unique_ptr<SimpleConfig> _config;
    std::unique_ptr<SimpleConfig> _bfConfig;

    const fhicl::ParameterSet       _simulatedDetector;

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
    bool _standardMu2eDetector;

    // All of the detectors that we know about.
    DetMap _detectors;

    // Keep a count of how many runs we have seen.
    int _run_count = 0;

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
