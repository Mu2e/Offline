#ifndef ConditionsService_ConditionsService_hh
#define ConditionsService_ConditionsService_hh

//
// Primitive conditions data service.
// It does not yet do validty checking.
//
// $Id: ConditionsService.hh,v 1.16 2012/07/15 22:06:16 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/07/15 22:06:16 $
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

// Framework include files
#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Services/Registry/ActivityRegistry.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "cetlib/exception.h"

// Mu2e include files.
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eInterfaces/inc/ConditionsEntity.hh"

// Other external include files.
#include "boost/shared_ptr.hpp"

namespace mu2e {

  class ConditionsService {

  public:
    ConditionsService(const fhicl::ParameterSet&, art::ActivityRegistry&);
    ~ConditionsService();

    void preBeginRun(art::Run const &);

    // Not sure if we really want this.  It might be abused more than used?
    SimpleConfig const& config() const { return _config;}

  private:

    // The name of the input file.  Later will be a db key or similar.
    std::string _conditionsFile;

    // For how the conditions data is held in the file managed by
    // this config object.  It can later evolve to a database.
    SimpleConfig _config;

    // Perform any consistency checks.
    void checkConsistency();

    typedef boost::shared_ptr<ConditionsEntity> ConditionsEntityPtr;
    typedef std::map<std::string,ConditionsEntityPtr> ConditionsMap;

    template <typename ENTITY> friend class ConditionsHandle;

    // For now the key and version arguments are ignored.
    template <class ENTITY>
    ENTITY* getElement( std::string const& , std::string const& )
    {
      if(_run_count==0)
        throw cet::exception("GEOM")
          << "Cannot get _entities from an unconfigured conditions service.\n"
          << "You've attempted to a get an element before the first run\n";

      // to use this generic way requires a map of names (typeid?) to
      // abstract elements.
      // find the conditions entity element requested
      std::string name = typeid(ENTITY).name();
      ConditionsMap::iterator it(_entities.find(name));
      if(it==_entities.end())
        throw cet::exception("GEOM")
          << "Failed to retrieve conditions entity of type " << name << "\n";

      // this must succeed or there is something terribly wrong
      ENTITY* d = dynamic_cast<ENTITY*>(it->second.get());

      if(d==0)
        throw cet::exception("GEOM")
          << "Failed to convert found conditions entity " << name
          << " to its correct type.  There is a serious problem.\n";

      return d;
    }


    ConditionsMap _entities;
    int _run_count;

    // This is not copyable or assignable - private and unimplemented.
    ConditionsService const& operator=(ConditionsService const& rhs);
    ConditionsService(ConditionsService const& rhs);

    // Don't need to expose definition of private template in header
    template <typename ENTITY> void addEntity(std::auto_ptr<ENTITY> d);

  };

}

#endif /* ConditionsService_ConditionsService_hh */
