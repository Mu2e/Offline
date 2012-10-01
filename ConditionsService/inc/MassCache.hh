#ifndef ConditionsService_MassCache_hh
#define ConditionsService_MassCache_hh
//
// Cache values of masses by PDG Id to reduce expensive
// lookups in the particle data table.
//
// $Id: MassCache.hh,v 1.2 2012/10/01 18:19:29 genser Exp $
// $Author: genser $
// $Date: 2012/10/01 18:19:29 $
//
// Original author Rob Kutschke
//

#include "MCDataProducts/inc/PDGCode.hh"

// C++ includes.
#include <iostream>
#include <map>

namespace mu2e {

  class MassCache {

  public:

    MassCache ();
    // Accept compiler generated:
    // copy c'tor, d'tor and assignment operator.

    typedef PDGCode::type id_type;

    double mass( id_type pdgId );
    size_t size() const { return cache_.size(); }

  private:

    typedef std::map<id_type,double>    map_type;
    typedef map_type::value_type      value_type;


    map_type cache_;
    double   lastMass_;
    id_type  lastId_;

    // Helper functions./
    double getMassFromPDT ( id_type id );

  };

  // Shift left (printing) operator.
  inline std::ostream& operator<<(std::ostream& ost,
                                  const MassCache& daqpar ){

    ost << "[ ";
    /*
    for ( map_type::const_iterator i=cache_.begin(), e=cache_.end();
          i !=e; ++i ){
      ost << "( "
          << i->first << ", "
          << i->second
          << " ) ";
    }
    */
    ost << " ]";
    return ost;
  }
}

#endif /* ConditionsService_MassCache_hh */
