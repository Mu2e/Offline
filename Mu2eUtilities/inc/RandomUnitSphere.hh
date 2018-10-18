#ifndef Mu2eUtilities_RandomUnitSphere_hh
#define Mu2eUtilities_RandomUnitSphere_hh

//
// Return CLHEP::Hep3Vector objects that are unit vectors uniformly
// distributed over the unit sphere.
//
// $Id: RandomUnitSphere.hh,v 1.11 2013/03/21 20:33:55 gandr Exp $
// $Author: gandr $
// $Date: 2013/03/21 20:33:55 $
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

  struct RandomUnitSphereParams {
    double czmin;
    double czmax;
    double phimin;
    double phimax;
    RandomUnitSphereParams(double ci, double ca, double pi, double pa)
      : czmin(ci), czmax(ca), phimin(pi), phimax(pa) {}
  };

  class RandomUnitSphere {
  public:

    explicit RandomUnitSphere( CLHEP::HepRandomEngine& engine,
                               double czmin=-1.,
                               double czmax=1.,
                               double phimin=0.,
                               double phimax=CLHEP::twopi);

    explicit RandomUnitSphere( CLHEP::HepRandomEngine& engine, const RandomUnitSphereParams& par);

    // I have inspected RandFlat source and assert that its destructor does not throw
    ~RandomUnitSphere() noexcept = default;

    CLHEP::Hep3Vector fire();

    // Alternate fire syntax which modifies the magnitude of the vector.
    CLHEP::Hep3Vector fire( double magnitude ){
      return magnitude*fire();
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

    CLHEP::HepRandomEngine& engine(){
      return _randFlat.engine();
    }

  private:

    double _czmin;
    double _czmax;
    double _phimin;
    double _phimax;

    // The underlying uniform random number distribution.
    CLHEP::RandFlat _randFlat;
  };

}

#endif /* Mu2eUtilities_RandomUnitSphere_hh */
