// Geometry of the hall, dirt, etc.
//
// Andrei Gaponenko, 2011

#ifndef MU2EHALL_HH
#define MU2EHALL_HH

#include <map>
#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/TwoVector.h"

#include "GeomPrimitives/inc/ExtrudedSolid.hh"
#include "Mu2eInterfaces/inc/Detector.hh"

#include "art/Persistency/Common/Wrapper.h"

namespace mu2e {

  class Mu2eHallMaker;
  
  class Mu2eHall : virtual public Detector {
  public:
    
    const std::map<std::string,ExtrudedSolid>& getBldgSolids() const { return bldgSolids_; }
    const std::map<std::string,ExtrudedSolid>& getDirtSolids() const { return dirtSolids_; }
    
    const ExtrudedSolid& 
    getBldgSolid( const std::string& str ) const {
      return bldgSolids_.find( str )->second;
    }

    const ExtrudedSolid& 
    getDirtSolid( const std::string& str ) const {
      return dirtSolids_.find( str )->second;
    }
 

    //----------------------------------------------------------------
  private:
    friend class Mu2eHallMaker;

    // Private ctr: the class should be only obtained via the maker
    Mu2eHall(){}

    // Or read back from persistent storage
    template<class T> friend class art::Wrapper;

    std::map<std::string,ExtrudedSolid> bldgSolids_;
    std::map<std::string,ExtrudedSolid> dirtSolids_;

  };

}

#endif/*MU2EHALL_HH*/
