#ifndef RANDOMUNITSHPERE_HH
#define RANDOMUNITSHPERE_HH

//
// Return CLHEP::Hep3Vector objects that are unit vectors uniformly
// distributed over the unit sphere.
// 
// $Id: RandomUnitSphere.hh,v 1.6 2010/08/22 20:09:52 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/08/22 20:09:52 $
//
// Original author Rob Kutschke
//
// Allows for range limitations on cos(theta) and phi.
//
//

#include <cmath>

#include "CLHEP/Random/RandFlat.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Units/PhysicalConstants.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e { 

  class RandomUnitSphere {

  public:

    explicit RandomUnitSphere( double czmin=-1.,
                               double czmax=1.,
                               double phimin=0.,
                               double phimax=CLHEP::twopi);

    explicit RandomUnitSphere( CLHEP::HepRandomEngine& engine,
                               double czmin=-1.,
                               double czmax=1.,
                               double phimin=0.,
                               double phimax=CLHEP::twopi);

    ~RandomUnitSphere(){}

    CLHEP::Hep3Vector fire() const;
    CLHEP::Hep3Vector fire( double magnitude ) const;

    // Depracated
    CLHEP::Hep3Vector shoot() const{
      return fire();
    }


    void setczmin(double czmin){
      _czmin=czmin;
    }

    void setczmax(double czmax){
      _czmax=czmax;
    }

    void setphimin(double phimin){
      _phimin=phimin;
    }

    void setphimax(double phimax){
      _phimax=phimax;
    }

    double czmin(){ return _czmin;}
    double czmax(){ return _czmax;}

    double phimin(){ return _phimin;}
    double phimax(){ return _phimax;}

  private:

    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    // The underlying uniform random number distribution.
    mutable CLHEP::RandFlat _randFlat;

  };

}

#endif
