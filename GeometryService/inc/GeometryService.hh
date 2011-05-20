#ifndef GeometryService_GeometryService_hh
#define GeometryService_GeometryService_hh

//
// Maintain up to date geometry information and serve it to
// other services and to the modules.
//
// $Id: GeometryService.hh,v 1.9 2011/05/20 12:23:42 greenc Exp $
// $Author: greenc $
// $Date: 2011/05/20 12:23:42 $
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

#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeometryService/inc/Detector.hh"
#include "boost/shared_ptr.hpp"

namespace mu2e {

// Forward declarations
  class Target;

  class GeometryService {
public:
    GeometryService(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~GeometryService();

    void preBeginRun( art::Run const &run);

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

    std::auto_ptr<SimpleConfig> _config;

    // Check the configuration.
    void checkConfig();

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

    // The name of the input file.  Later will be a db key or similar.
    std::string _inputfile;

    DetMap _detectors;
    int _run_count;

    // This is not copyable or assignable - private and unimplemented.
    GeometryService const& operator=(GeometryService const& rhs);
    GeometryService(GeometryService const& rhs);

    template <typename DET>
    void addDetector(std::auto_ptr<DET> d)
    {
      if(_detectors.find(typeid(DET).name())!=_detectors.end())
        throw cet::exception("GEOM") << "failed to install detector "
                                     << d->name() << "\nwith type name "
                                     << typeid(DET).name() << "\n";

      DetectorPtr ptr(d.release());
      _detectors[typeid(DET).name()] = ptr;
    }


  };

}

#endif /* GeometryService_GeometryService_hh */
