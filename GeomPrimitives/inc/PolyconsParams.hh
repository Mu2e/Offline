#ifndef GeomPrimitives_PolyconsParams_hh
#define GeomPrimitives_PolyconsParams_hh

//
// The parameters of a Polycone
//
// $Id: PolyconsParams.hh,v 1.1 2013/04/30 14:56:57 knoepfel Exp $
// $Author: knoepfel $
// $Date: 2013/04/30 14:56:57 $
//
// Original author KLG
//

#include <vector>
#include <string>
#include <ostream>

#include "CLHEP/Units/PhysicalConstants.h"

namespace mu2e {

  class PolyconsParams{

  public:

    PolyconsParams(const std::vector<double>& zPlanes,
                   const std::vector<double>& rInner,
                   const std::vector<double>& rOuter,
                   double phiStart = 0.,
                   double phiTotal = CLHEP::twopi) :
      _zPlanes( zPlanes ), _rInner( rInner ), _rOuter( rOuter ),
      _phiStart( phiStart ), _phiTotal( phiTotal )
    {}

    // Use compiler-generated copy c'tor, copy assignment, and d'tor.

    double phi0()           const { return _phiStart; }
    double phiTotal()       const { return _phiTotal; }

    const std::vector<double>& zPlanes() const { return _zPlanes; }
    const std::vector<double>& rInner()  const { return _rInner; }
    const std::vector<double>& rOuter()  const { return _rOuter; }

    unsigned numZPlanes() const { return _zPlanes.size(); }

  private:

    std::vector<double> _zPlanes;
    std::vector<double> _rInner;
    std::vector<double> _rOuter;

    double _phiStart;
    double _phiTotal;

  };

  inline std::ostream& operator<<(std::ostream& ost,
                                  const PolyconsParams& pp ){
    ost << "(";
    for (  std::vector<double>::size_type i=0; i !=pp.zPlanes().size() ; ++i ) {
        ost << pp.zPlanes()[i] << ", "
            << pp.rInner()[i]  << ", "
            << pp.rOuter()[i]  << "; ";
    }
    ost << pp.phi0()   << " "
        << pp.phiTotal() << " )";
    return ost;
  }

}

#endif /* GeomPrimitives_PolyconsParams_hh */
