#ifndef BFieldGeom_BFMapBase_hh
#define BFieldGeom_BFMapBase_hh
//
// Interface to the magnetic field maps. Used by BFMap and BFMapSet.
//
// $Id: BFMapBase.hh,v 1.6 2011/05/18 02:27:14 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:14 $
//

#include <string>
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class BFMapType;

  class BFMapBase {

  public:

    virtual ~BFMapBase();

    virtual bool getBFieldWithStatus(const CLHEP::Hep3Vector &,
				     CLHEP::Hep3Vector & ) const = 0;

    CLHEP::Hep3Vector getBField(const CLHEP::Hep3Vector& pos) const {
      CLHEP::Hep3Vector result;
      if( getBFieldWithStatus(pos,result) ) return result;
      else return CLHEP::Hep3Vector(0,0,0);
    }

    virtual const std::string& getKey() const = 0;
    virtual bool isValid(CLHEP::Hep3Vector const& point) const = 0;

    virtual BFMapType type() const = 0;

  };

} // end namespace mu2e

#endif /* BFieldGeom_BFMapBase_hh */
