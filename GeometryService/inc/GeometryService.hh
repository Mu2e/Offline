#ifndef GeometryService_GeometryService_hh
#define GeometryService_GeometryService_hh

//
// Maintain up to date geometry information and serve it to
// other services and to the modules.
//
// $Id: GeometryService.hh,v 1.2 2010/02/11 15:45:43 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/02/11 15:45:43 $
//
// Original author Rob Kutschke
//

// C++ include files
#include <string>
#include <memory>

// Framework include files
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/ActivityRegistry.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "GeometryService/inc/Detector.hh"
#include "boost/shared_ptr.hpp"

namespace mu2e {

// Forward declarations
  class Target;
  class CTracker;

  class GeometryService {
public:
    GeometryService(const edm::ParameterSet&, edm::ActivityRegistry&);
    ~GeometryService();
    
    void preBeginRun( edm::RunID const& id, edm::Timestamp const& ts);

    SimpleConfig const& config() const { return *_config;}

    template <class DET>
    bool hasElement()
    {
      if(_run_count==0) 
	throw cms::Exception("GEOM")
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
	throw cms::Exception("GEOM")
	  << "Cannot get detectors from an unconfigured geometry service.\n"
	  << "You've attempted to a get an element before the first run\n";
      
      // to use this generic way requires a map of names (typeid?) to
      // abstract elements.
      // find the detector element requested
      std::string name = typeid(DET).name();
      DetMap::iterator it(_detectors.find(name));
      if(it==_detectors.end())
	throw cms::Exception("GEOM")
	  << "Failed to retrieve detector element of type " << name << "\n";

      // this must succeed or there is something terribly wrong
      DET* d = dynamic_cast<DET*>(it->second.get());

      if(d==0)
	cms::Exception("GEOM")
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
	throw cms::Exception("GEOM") << "failed to install detector "
				     << d->name() << "\nwith type name "
				     << typeid(DET).name() << "\n";
      
      DetectorPtr ptr(d.release());
      _detectors[typeid(DET).name()] = ptr;
    }
    

  };

}

#endif
