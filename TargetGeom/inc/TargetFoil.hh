#ifndef TargetGeom_TargetFoil_hh
#define TargetGeom_TargetFoil_hh

//
// Class to represent one target foil.
// For now these are just disks perpendicular to the z axis.
//
// $Id: TargetFoil.hh,v 1.7 2012/08/27 22:19:43 mf Exp $
// $Author: mf $
// $Date: 2012/08/27 22:19:43 $
//
// Original author Rob Kutschke
//
// Coordinates are given in the detector coordinate
// system in cm.
//

#include <string>

// Includes from CLHEP
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TargetFoil {

  public:
    TargetFoil( int id,
                CLHEP::Hep3Vector const& c,
                CLHEP::Hep3Vector const& n,
                double rOut,
                double rIn,
                double t,
                std::string m):
      _id(id),
      _c(c),
      _norm(n),
      _rOut(rOut),
      _rIn(rIn),
      _t(t),
      _material(m){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    int id() const { return _id; }

    CLHEP::Hep3Vector const& center()  const { return _c;}
    CLHEP::Hep3Vector const& normal()  const { return _norm;}

    double rOut()          const { return _rOut;}
    double rIn()           const { return _rIn;}
    double halfThickness() const { return _t;}

    std::string material() const { return _material;}

  private:

    int _id;

    // Center of the foil.
    CLHEP::Hep3Vector _c;
    // "+z" normal vector
    CLHEP::Hep3Vector _norm;

    // Inner and outer radii.
    double _rOut;
    double _rIn;

    // Thickness in z.
    double _t;
                // Is this the thickness or half thickness??
                
    std::string _material;

  };

}
#endif /* TargetGeom_TargetFoil_hh */
