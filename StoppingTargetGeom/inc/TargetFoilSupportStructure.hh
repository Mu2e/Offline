#ifndef StoppingTargetGeom_TargetFoilSupportStructure_hh
#define StoppingTargetGeom_TargetFoilSupportStructure_hh

//
// Class to represent one support structure of one target foil.
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

  class TargetFoilSupportStructure {

  public:
    TargetFoilSupportStructure( int support_id, int foil_id,
                CLHEP::Hep3Vector const& foilCenterInMu2e,
                CLHEP::Hep3Vector const& n,
                double radius,
                double length,
                double angleOffset,
                double foil_outer_radius,
                std::string const& m,
                const CLHEP::Hep3Vector& detSysOrigin
                ):
      _support_id(support_id),
      _foil_id(foil_id),
      _centerInMu2e(foilCenterInMu2e),
      _centerInDetectorSystem(foilCenterInMu2e -  detSysOrigin),
      _norm(n),
      _radius(radius),
      _length(length),
      _angleOffset(angleOffset),
      _foil_outer_radius(foil_outer_radius),
      _material(m){
    }

    // Use compiler-generated copy c'tor, copy assignment, and d'tor

    int support_id() const { return _support_id; }
    int foil_id() const { return _foil_id; }

    CLHEP::Hep3Vector const& centerInMu2e()  const { return _centerInMu2e;}
    CLHEP::Hep3Vector const& centerInDetectorSystem()  const { return _centerInDetectorSystem;}

    CLHEP::Hep3Vector const& normal()  const { return _norm;}

    double radius()          const { return _radius;}
    double length()          const { return _length;}
    double angleOffset()          const { return _angleOffset;}
    double foil_outer_radius()          const { return _foil_outer_radius;}

    std::string material() const { return _material;}

  private:

    int _support_id;
    int _foil_id;

    // Center of the foil.
    CLHEP::Hep3Vector _centerInMu2e;
    CLHEP::Hep3Vector _centerInDetectorSystem;

    // "+z" normal vector
    CLHEP::Hep3Vector _norm;

    double _radius; // radius of the supporting structure wire

    double _length; // length of the supporting structure wire
    double _angleOffset; // angle of first support wire wrt vertical
    double _foil_outer_radius; // outer radius of the foil the wire connects to

    std::string _material;

  };

}
#endif /* StoppingTargetGeom_TargetFoilSupportStructure_hh */
