#ifndef GeometryService_GeometryService_hh
#define GeometryService_GeometryService_hh

//
// Maintain up to date geometry information and serve it to
// other services and to the modules.
//
// $Id: GeometryService.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
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

    // A test method.
    std::string getInfo() const{ return "Here is a string"; }

#if 0
    // I don't like these access methods but they are what we have for now.
    bool hasTarget() const { return (target_ !=0);}
    bool hasCTracker() const { return (ctracker_ !=0);}
    // Some detector elements we have not yet defined.
    bool hasECal()       const { return false;}
    bool hasCosmicVeto() const { return false;}
    bool hasBeamLine()   const { return false;}
#endif

    SimpleConfig const& config() const { return *_config;}

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
