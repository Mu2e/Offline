#ifndef GeomPrimitives_PolyhedraParams_hh
#define GeomPrimitives_PolyhedraParams_hh

//
// The parameters of a Polyhedron
//
// Blatantly stolen from Kyle's PolyconsParams class.
// David Norvil Brown (UofL), September 2017
//
// Original author KLG
//

#include <vector>
#include <string>
#include <ostream>

#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  class PolyhedraParams{

  public:

    PolyhedraParams( int  nSides,
		     const std::vector<double>& zPlanes,
		     const std::vector<double>& rInner,
		     const std::vector<double>& rOuter,
		     double phiStart = 0.,
		     double phiTotal = CLHEP::twopi) :
      _nSides( nSides ), _zPlanes( zPlanes ), _rInner( rInner ), 
      _rOuter( rOuter ),
      _phiStart( phiStart ), _phiTotal( phiTotal )
    {}

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    double phi0()           const { return _phiStart; }
    double phiTotal()       const { return _phiTotal; }
    double nSides()         const { return _nSides;   }

    const std::vector<double>& zPlanes() const { return _zPlanes; }
    const std::vector<double>& rInner()  const { return _rInner; }
    const std::vector<double>& rOuter()  const { return _rOuter; }

    unsigned numZPlanes() const { return _zPlanes.size(); }

  private:

    int                 _nSides;

    std::vector<double> _zPlanes;
    std::vector<double> _rInner;
    std::vector<double> _rOuter;

    double              _phiStart;
    double              _phiTotal;

  };

}

#endif /* GeomPrimitives_PolyhedraParams_hh */
