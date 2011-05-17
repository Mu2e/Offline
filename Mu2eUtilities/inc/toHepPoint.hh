#ifndef Mu2eUtilities_toHepPoint_hh
#define Mu2eUtilities_toHepPoint_hh
//
//  Free functions to convert between  HepPoint and CLHEP::Hep3Vector.
//
// $Id: toHepPoint.hh,v 1.3 2011/05/17 15:41:36 greenc Exp $
// $Author: greenc $ 
// $Date: 2011/05/17 15:41:36 $
//
// Original author Rob Kutschke
//

// Modern CLHEP
#include "CLHEP/Vector/ThreeVector.h"

// Ancient CLHEP, copied from BaBar
#include "CLHEP/Geometry/HepPoint.h"

namespace mu2e {

  inline const HepPoint toHepPoint( const CLHEP::Hep3Vector& v){
    return HepPoint( v.x(), v.y(), v.z());
  }

  inline const CLHEP::Hep3Vector fromHepPoint( const HepPoint& v){
    return CLHEP::Hep3Vector( v.x(), v.y(), v.z() );
  }

}
#endif /* Mu2eUtilities_toHepPoint_hh */
