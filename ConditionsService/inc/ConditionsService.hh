#ifndef ConditionsService_HH
#define ConditionsService_HH

//
// Primitive conditions data service.
// It does not yet do validty checking.
//
// $Id: ConditionsService.hh,v 1.2 2009/11/12 01:35:50 kutschke Exp $
// $Author: kutschke $
// $Date: 2009/11/12 01:35:50 $
//
// Original author Rob Kutschke
//
//
// Notes
// 1) There are two types of accessors.
//     a) ones that just forward to the config file, for those types of entities.
//     b) The method that we want to move to eventually.
//        CalibHandle<T> handle;
//        conditionsService->get( key, handle, version)      
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
#include "ConditionsService/inc/ConditionsEntity.hh"
#include "boost/shared_ptr.hpp"

namespace mu2e {

  // Forward declarations

  class ConditionsService {

  public:
    ConditionsService(const edm::ParameterSet&, edm::ActivityRegistry&);
    ~ConditionsService();
    
    void preBeginRun( edm::RunID const& id, edm::Timestamp const& ts);
    

    // There are two types of accessors.
    // For
    SimpleConfig const& config() const { return *_config;}

    bool getBool ( std::string const& key){
      return _config->getBool(key);
    }

    int getInt( std::string const& key){
      return _config->getInt(key);
    }

    double getDouble( std::string const& key){
      return _config->getDouble(key);
    }

  private:
    
    // For how the conditions data is held in the file managed by
    // this config object.  It can later evolve to a database.
    std::auto_ptr<SimpleConfig> _config;

    // Perform any consistency checks.
    void checkConsistency();
    
    typedef boost::shared_ptr<ConditionsEntity> ConditionsEntityPtr;
    typedef std::map<std::string,ConditionsEntityPtr> ConditionsMap;
    
    template <typename ENTITY> friend class ConditionsHandle;
    
    // For now the key and version arguments are ignored.
    template <class ENTITY>
    ENTITY* getElement( std::string const& key, std::string const& version)
    {
      if(_run_count==0) 
	throw cms::Exception("GEOM")
	  << "Cannot get _entities from an unconfigured conditions service.\n"
	  << "You've attempted to a get an element before the first run\n";
      
      // to use this generic way requires a map of names (typeid?) to
      // abstract elements.
      // find the conditions entity element requested
      std::string name = typeid(ENTITY).name();
      ConditionsMap::iterator it(_entities.find(name));
      if(it==_entities.end())
	throw cms::Exception("GEOM")
	  << "Failed to retrieve conditions entity of type " << name << "\n";

      // this must succeed or there is something terribly wrong
      ENTITY* d = dynamic_cast<ENTITY*>(it->second.get());
      
      if(d==0)
	cms::Exception("GEOM")
	  << "Failed to convert found conditions entity " << name
	  << " to its correct type.  There is a serious problem.\n";
      
      return d;
    }
    
    // The name of the input file.  Later will be a db key or similar.
    std::string _conditionsFile;

    ConditionsMap _entities;
    int _run_count;

    // This is not copyable or assignable - private and unimplemented.
    ConditionsService const& operator=(ConditionsService const& rhs);
    ConditionsService(ConditionsService const& rhs);
    
    template <typename ENTITY>
    void addEntity(std::auto_ptr<ENTITY> d)
    {
      if(_entities.find(typeid(ENTITY).name())!=_entities.end())
	throw cms::Exception("GEOM") << "failed to install conditions entity "
				     << d->name() << "\nwith type name "
				     << typeid(ENTITY).name() << "\n";
      
      ConditionsEntityPtr ptr(d.release());
      _entities[typeid(ENTITY).name()] = ptr;
    }

  };

}

#endif
