#ifndef hep3VectorFromStdVector_HH
#define hep3VectorFromStdVector_HH

//
//  A variety of ways of copying a std::vector into a 
//  Hep3Vector.
//
// $Id: hep3VectorFromStdVector.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
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
#endif
