#ifndef SeedService_EngineId_hh
#define SeedService_EngineId_hh
//
// An identifier for random engines.  An identifier may consist
// of simply a module label or a module label plus an instance name.
//
//
// Contact person Rob Kutschke
//

#include <string>
#include <tuple>
#include <iostream>

namespace mu2e {

  namespace SeedServiceHelper {

    struct EngineId{

      EngineId( std::string const& mod, std::string const& inst):
        moduleLabel(mod),
        instanceName(inst),
        instanceDefined(true){}

      EngineId( std::string const& mod):
        moduleLabel(mod),
        instanceName(),
        instanceDefined(false){}

      // Accept compiler written d'tor, copy c'tor and copy assignment.

      // An id without an instance name is distinct from every id with one,
      // so that each engine of a module receives its own seed.
      bool operator==( EngineId const& rhs ) const{
        return std::tie(moduleLabel, instanceDefined, instanceName) ==
          std::tie(rhs.moduleLabel, rhs.instanceDefined, rhs.instanceName);
      }

      bool operator<( EngineId const& rhs ) const{
        return std::tie(moduleLabel, instanceDefined, instanceName) <
          std::tie(rhs.moduleLabel, rhs.instanceDefined, rhs.instanceName);
      }

      std::string moduleLabel;
      std::string instanceName;
      bool instanceDefined;

    }; // end class EngineId

    inline std::ostream& operator<<(std::ostream& ost,
                                    const EngineId& id ){
      ost << id.moduleLabel;
      if ( id.instanceDefined ){
        ost << "." << id.instanceName;
      }
      return ost;
    }

  } // end namespace SeedServiceHelper

} // end namespace mu2e

#endif /* SeedService_EngineId_hh */
