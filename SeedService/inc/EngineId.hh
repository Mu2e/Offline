#ifndef SeedService_EngineId_hh
#define SeedService_EngineId_hh
//
// An identifier for random engines.  An identifier may consist
// of simply a module label or a module label plus an instance name.
//
// $Id: EngineId.hh,v 1.1 2012/02/04 00:12:48 kutschke Exp $
// $Author: kutschke $
// $Date: 2012/02/04 00:12:48 $
//
// Contact person Rob Kutschke
//

#include <string>
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

      bool operator==( EngineId const& rhs ) const{
        if ( moduleLabel  != rhs.moduleLabel  ) return false;
        if ( instanceDefined && rhs.instanceDefined ) {
          if ( instanceName != rhs.instanceName ) return false;
        }
        return true;
      }

      bool operator<( EngineId const& rhs ) const{
        if ( moduleLabel  < rhs.moduleLabel  ) return true;
        if ( instanceDefined && rhs.instanceDefined ) {
          if ( moduleLabel == rhs.moduleLabel  ) {
            if ( instanceName < rhs.instanceName ) return true;
          }
        }
        return false;
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
