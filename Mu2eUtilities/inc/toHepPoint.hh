#ifndef Mu2eUtilities_toHepPoint_HH
#define Mu2eUtilities_toHepPoint_HH
//
//  Free functions to convert between  HepPoint and CLHEP::Hep3Vector.
//
// $Id: toHepPoint.hh,v 1.2 2010/05/17 21:47:32 genser Exp $
// $Author: genser $ 
// $Date: 2010/05/17 21:47:32 $
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
#endif
