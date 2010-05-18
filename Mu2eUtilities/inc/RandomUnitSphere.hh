#ifndef RANDOMUNITSHPERE_HH
#define RANDOMUNITSHPERE_HH

//
// Return CLHEP::Hep3Vector objects that are unit vectors uniformly
// distributed over the unit sphere.
// 
// $Id: RandomUnitSphere.hh,v 1.3 2010/05/18 20:28:50 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2010/05/18 20:28:50 $
//
// Original author Rob Kutschke
//
// Allows for range limitations on cos(theta) and phi.
//
//

#include <cmath>
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e { 

  class RandomUnitSphere {

  public:
  
    RandomUnitSphere():
      _czmin(-1.),
      _czmax( 1.),
      _phimin(0.),
      _phimax(CLHEP::twopi){
    }

    RandomUnitSphere( double czmin,
                      double czmax,
                      double phimin=0,
                      double phimax=CLHEP::twopi):
      _czmin(czmin),
      _czmax(czmax),
      _phimin(phimin),
      _phimax(phimax){
    }

    ~RandomUnitSphere(){}

    CLHEP::Hep3Vector shoot() const;

    
    // Do I really want the setters?
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
    
  };

}

#endif
