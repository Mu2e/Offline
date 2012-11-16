#ifndef GeometryService_GeometryService_hh
#define GeometryService_GeometryService_hh

//
// Maintain up to date geometry information and serve it to
// other services and to the modules.
//
// $Id: GeometryService.hh,v 1.19 2012/11/16 23:30:51 genser Exp $
// $Author: genser $
// $Date: 2012/11/16 23:30:51 $
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
#include "cetlib/exception.h"

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

  class GeometryService {
public:
    GeometryService(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~GeometryService();

    // Functions registered for callbacks.
    void preBeginRun( art::Run const &run);
    void postEndJob();

    SimpleConfig const& config() const { return *_config;}

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
    std::auto_ptr<SimpleConfig> _config;

    // Check the configuration.
    void checkConfig();
    void checkTrackerConfig();
    void checkCalorimeterConfig();

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

    // All of the detectors that we know about.
    DetMap _detectors;

    // Keep a count of how many runs we have seen.
    int _run_count;

    // This is not copyable or assignable - private and unimplemented.
    GeometryService const& operator=(GeometryService const& rhs);
    GeometryService(GeometryService const& rhs);

    // Don't need to expose definition of private template in header
    template <typename DET> void addDetector(std::auto_ptr<DET> d);
    template <typename DETALIAS, typename DET> void addDetectorAliasToBaseClass(std::auto_ptr<DET> d);

    // Some information that is provided through the GeometryService
    // should only be used inside GEANT jobs.  The following method is
    // used by G4 to make this info available.
    friend class G4;
    friend class Mu2eG4Study;
    void addWorldG4();

  };

}

#endif /* GeometryService_GeometryService_hh */
