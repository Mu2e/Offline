#ifndef GeomPrimitives_Polyhedra_hh
#define GeomPrimitives_Polyhedra_hh

//
// The parameters of a Polyhedra
//
// This is a blatant rip-off of Kyle's Polycone class, but now
// used for representing polyhedra, which correlates with 
// Geant4's G4Polyhedra class, without the specific dependence on G4.
// David Norvil Brown, September 2017
//
// Original author KLG
//

#include <vector>
#include <string>
#include <ostream>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "GeomPrimitives/inc/PolyhedraParams.hh"

namespace mu2e {

  class Polyhedra{

    //   G4Polyhedra( const G4String& name,
    //               G4double phiStart,     // initial phi starting angle
    //               G4double phiTotal,     // total phi angle
    //               G4int numZPlanes,     // number of z planes
    //               const G4double zPlane[],  // position of z planes
    //               const G4double rInner[],  // tangent distance to inner surface
    //               const G4double rOuter[])  // tangent distance to outer surface

  public:

    Polyhedra( int nSides,
	       const std::vector<double>& zPlanes,
	       const std::vector<double>& rInner,
	       const std::vector<double>& rOuter,
	       const CLHEP::Hep3Vector& originInMu2e,
	       const std::string& materialName,
	       double phiStart = 0.,
	       double phiTotal = CLHEP::twopi);

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    double phi0()           const { return _params.phi0(); }
    double phiTotal()       const { return _params.phiTotal(); }
    int    nSides()         const { return _params.nSides(); }

    const std::string& materialName() const { return _materialName; }

    const CLHEP::Hep3Vector& originInMu2e() const { return _originInMu2e; }

    const std::vector<double>& zPlanes() const { return _params.zPlanes(); }
    const std::vector<double>& rInner()  const { return _params.rInner(); }
    const std::vector<double>& rOuter()  const { return _params.rOuter(); }

    unsigned numZPlanes() const { return _params.numZPlanes(); }

    const PolyhedraParams & getPolyhedraParams() const { return _params; }

  private:

    PolyhedraParams      _params;

    CLHEP::Hep3Vector _originInMu2e;
    std::string       _materialName;

  };

  std::ostream& operator<<(std::ostream& os, const Polyhedra& p);
}

#endif /* GeomPrimitives_Polyhedra_hh */
