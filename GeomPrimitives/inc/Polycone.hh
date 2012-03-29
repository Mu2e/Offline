#ifndef GeomPrimitives_Polycone_hh
#define GeomPrimitives_Polycone_hh

//
// The parameters of a Polycone
//
// $Id: Polycone.hh,v 1.3 2012/03/29 19:07:23 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/29 19:07:23 $
//
// Original author KLG
//

#include <vector>
#include <string>
#include <cmath>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  class Polycone{

    //   G4Polycone( const G4String& name,
    //               G4double phiStart,     // initial phi starting angle
    //               G4double phiTotal,     // total phi angle
    //               G4int numZPlanes,     // number of z planes
    //               const G4double zPlane[],  // position of z planes
    //               const G4double rInner[],  // tangent distance to inner surface
    //               const G4double rOuter[])  // tangent distance to outer surface

  public:

    Polycone(const std::vector<double>& zPlanes,
             const std::vector<double>& rInner,
             const std::vector<double>& rOuter,
             const CLHEP::Hep3Vector& originInMu2e,
             const std::string& materialName,
             double phiStart = 0.,
             double phiTotal = CLHEP::twopi);

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    double phi0()           const { return _phiStart; }
    double phiTotal()       const { return _phiTotal; }

    const std::string& materialName() const { return _materialName; }

    const CLHEP::Hep3Vector& originInMu2e() const { return _originInMu2e; }

    const std::vector<double>& zPlanes() const { return _zPlanes; }
    const std::vector<double>& rInner()  const { return _rInner; }
    const std::vector<double>& rOuter()  const { return _rOuter; }

    unsigned numZPlanes() const { return _zPlanes.size(); }

  private:

    std::vector<double> _zPlanes;
    std::vector<double> _rInner;
    std::vector<double> _rOuter;

    CLHEP::Hep3Vector _originInMu2e;
    std::string       _materialName;

    double _phiStart;
    double _phiTotal;

  // do we need/want the rotation here?

};

// TBD
//   inline std::ostream& operator<<(std::ostream& ost,
//                                   const Polycone& tp ){
//     ost << "("
//         << tp.innerRadius() << " "
//         << tp.outerRadius() << " "
//         << tp.zHalfLength() << " "
//         << tp.phi0()        << " "
//         << tp.phiMax()      << ")"
//         << " )";
//     return ost;
//   }

// could use something like
//      std::copy (rOuter[0],rOuter[_numZPlanes],
//                 ostream_iterator<double>(std::cout,", "));



}


#endif /* GeomPrimitives_TubsParams_hh */
