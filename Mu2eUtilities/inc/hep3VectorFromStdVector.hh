#ifndef Mu2eUtilities_hep3VectorFromStdVector_hh
#define Mu2eUtilities_hep3VectorFromStdVector_hh

//
//  A variety of ways of copying a std::vector into a
//  CLHEP::Hep3Vector.
//
// $Id: hep3VectorFromStdVector.hh,v 1.5 2011/05/18 02:27:18 wb Exp $
// $Author: wb $
// $Date: 2011/05/18 02:27:18 $
//
// Original author Rob Kutschke
//

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {


  // Setter style version
  inline void
  hep3VectorFromStdVector( CLHEP::Hep3Vector cv,
                           std::vector<double> const& sv){
    cv.set( sv[0], sv[1], sv[2] );
  }

  // return by value version.
  inline CLHEP::Hep3Vector
  hep3VectorFromStdVector( std::vector<double> const& sv){
    return CLHEP::Hep3Vector( sv[0], sv[1], sv[2] );
  }

  // Safe version of setter.
  inline void
  checkedHep3VectorFromStdVector( CLHEP::Hep3Vector cv,
                                  std::vector<double> const& sv ){
    cv.set( sv.at(0), sv.at(1), sv.at(2) );
  }

  // Safe version of return by value.
  inline CLHEP::Hep3Vector
  checkedHep3VectorFromStdVector( std::vector<double> const& sv ){
    return CLHEP::Hep3Vector( sv.at(0), sv.at(1), sv.at(2) );
  }


}
#endif /* Mu2eUtilities_hep3VectorFromStdVector_hh */
