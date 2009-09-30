#ifndef Target_TargetFoil_HH
#define Target_TargetFoil_HH

//
// Class to represent one target foil.
// For now these are just disks perpendicular to the z axis.
//
// $Id: TargetFoil.hh,v 1.1 2009/09/30 22:57:47 kutschke Exp $
// $Author: kutschke $ 
// $Date: 2009/09/30 22:57:47 $
//
// Original author Rob Kutschke
//
// Coordinates are given in the detector coordinate 
// system in cm.
//

// Includes from CLHEP
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TargetFoil {

  public:
    TargetFoil( int id,
		CLHEP::Hep3Vector const& c, 
		double rOut, 
		double rIn, 
		double t):
      _id(id),
      _c(c),
      _rOut(rOut),
      _rIn(rIn),
      _t(t){
    }
    ~TargetFoil(){};

    int id() const { return _id; }

    CLHEP::Hep3Vector const& center()  const { return _c;} 

    double rOut()          const { return _rOut;}
    double rIn()           const { return _rIn;}
    double halfThickness() const { return _t;}

  private:

    int _id;

    // Center of the foil.
    CLHEP::Hep3Vector _c;
    
    // Inner and outer radii.
    double _rOut;
    double _rIn;

    // Thickness in z.
    double _t;
    

};

}
#endif
