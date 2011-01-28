#ifndef BFMAPBASE_HH
#define BFMAPBASE_HH
//
// Interface to the magnetic field maps. Used by BFMap and BFMapSet.
//
// $Id: BFMapBase.hh,v 1.2 2011/01/28 23:51:59 kutschke Exp $
// $Author: kutschke $
// $Date: 2011/01/28 23:51:59 $
//

#include <string>
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class BFMapBase {

  public:
    virtual ~BFMapBase();

    virtual CLHEP::Hep3Vector getBField(CLHEP::Hep3Vector const&) const = 0;
    virtual const std::string& getKey() const = 0;
    virtual bool isValid(CLHEP::Hep3Vector const& point) const = 0;

  };

} // end namespace mu2e

#endif
