#ifndef Mu2eUtilities_toHepPoint_HH
#define Mu2eUtilities_toHepPoint_HH
//
//  Free functions to convert between  HepPoint and Hep3Vector.
//
// $Id: toHepPoint.hh,v 1.1 2010/04/23 20:07:29 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/04/23 20:07:29 $
//
// Original author Rob Kutschke
//

// Modern CLHEP
#include "CLHEP/Vector/ThreeVector.h"

// Ancient CLHEP, copied from BaBar
#include "CLHEP/Geometry/HepPoint.h"

namespace mu2e {

  inline const HepPoint toHepPoint( const Hep3Vector& v){
    return HepPoint( v.x(), v.y(), v.z());
  }

  inline const Hep3Vector fromHepPoint( const HepPoint& v){
    return Hep3Vector( v.x(), v.y(), v.z() );
  }

}
#endif
