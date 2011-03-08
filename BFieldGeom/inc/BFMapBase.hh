#ifndef BFMAPBASE_HH
#define BFMAPBASE_HH
//
// Interface to the magnetic field maps. Used by BFMap and BFMapSet.
//
// $Id: BFMapBase.hh,v 1.4 2011/03/08 00:40:23 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/03/08 00:40:23 $
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

#endif
