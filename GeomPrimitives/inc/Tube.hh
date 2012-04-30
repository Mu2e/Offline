#ifndef GeomPrimitives_Tube_hh
#define GeomPrimitives_Tube_hh
//
// "Typical" Tube object
//
// $Id: Tube.hh,v 1.4 2012/04/30 16:21:31 gandr Exp $
// $Author: gandr $
// $Date: 2012/04/30 16:21:31 $
//
// Original author KLG
//

#include <string>
#include "cpp0x/array"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

#include "GeomPrimitives/inc/TubsParams.hh"


namespace mu2e {

  class TubsParams;
  class Tube {

  public:

    Tube(double rIn,double rOut,double halfLength, double phi0, double phiMax,
         std::string const & materialName, CLHEP::Hep3Vector const & originInMu2e);

    Tube(std::string const & materialName, CLHEP::Hep3Vector const & originInMu2e,
         double rIn,double rOut,double halfLength, 
         double phi0 = 0., double phiMax = CLHEP::twopi);

    double innerRadius() const { return _params.innerRadius(); }
    double outerRadius() const { return _params.outerRadius(); }
    double halfLength()  const { return _params.zHalfLength(); }
    double zHalfLength() const { return _params.zHalfLength(); }
    double phi0()        const { return _params.phi0(); }
    double phiMax()      const { return _params.phiMax(); }

    std::string const & materialName() const { return _materialName; }

    CLHEP::Hep3Vector const & originInMu2e() const { return _originInMu2e; }

    TubsParams const & getTubsParams() const { return _params; }

    // Genreflex can't do persistency of vector<Tube> without a default constructor
    Tube() : _params(0,0,0) {}

  private:

    TubsParams        _params;
    CLHEP::Hep3Vector _originInMu2e;
    std::string       _materialName;

    // do we need/want the rotation here?

  };

}

#endif/*GeomPrimitives_Tube_hh*/
