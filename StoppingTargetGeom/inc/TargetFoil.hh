#ifndef StoppingTargetGeom_TargetFoil_hh
#define StoppingTargetGeom_TargetFoil_hh

//
// Class to represent one target foil.
// For now these are just disks perpendicular to the z axis.
//
//
// Original author Rob Kutschke
//
// Coordinates are given in the detector coordinate
// system in mm.
//

#include <string>

// Includes from CLHEP
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class TargetFoil {

  public:
    TargetFoil( int id,
                CLHEP::Hep3Vector const& foilCenterInMu2e,
                CLHEP::Hep3Vector const& n,
                double rOut,
                double rIn,
                double t,
                std::string const& m,
                const CLHEP::Hep3Vector& detSysOrigin
                ):
      _id(id),
      _centerInMu2e(foilCenterInMu2e),
      _centerInDetectorSystem(foilCenterInMu2e -  detSysOrigin),
      _norm(n),
      _rOut(rOut),
      _rIn(rIn),
      _t(t),
      _material(m){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    int id() const { return _id; }

    CLHEP::Hep3Vector const& centerInMu2e()  const { return _centerInMu2e;}
    CLHEP::Hep3Vector const& centerInDetectorSystem()  const { return _centerInDetectorSystem;}

    CLHEP::Hep3Vector const& normal()  const { return _norm;}

    double rOut()          const { return _rOut;}
    double rIn()           const { return _rIn;}
    double halfThickness() const { return _t;}

    std::string material() const { return _material;}

  private:

    int _id;

    // Center of the foil.
    CLHEP::Hep3Vector _centerInMu2e;
    CLHEP::Hep3Vector _centerInDetectorSystem;

    // "+z" normal vector
    CLHEP::Hep3Vector _norm;

    // Inner and outer radii.
    double _rOut;
    double _rIn;

    // Half-thickness in z.
    double _t;

    std::string _material;

  };

}
#endif /* StoppingTargetGeom_TargetFoil_hh */
