#ifndef CosmicRayShieldGeom_CRSAbsorberLayer_hh
#define CosmicRayShieldGeom_CRSAbsorberLayer_hh
//
// Representation of one Absorber Layer in  CosmicRayShield
//
//

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e 
{

  class CRSAbsorberLayer
  {

    friend class CRSScintillatorModule;
    friend class CRSScintillatorShield;
    friend class CosmicRayShieldMaker;

    public:

    CRSAbsorberLayer() {}

    CRSAbsorberLayer(const CLHEP::Hep3Vector &position, const std::vector<double> &halfLength);
    // Accept the compiler generated destructor, copy constructor and assignment operators

    const CLHEP::Hep3Vector &getPosition() const {return _position;}
    const std::vector<double> &getHalfLengths() const {return _halfLengths;}

    private:

    CLHEP::Hep3Vector _position;
    std::vector<double> _halfLengths;
  };
}

#endif /* CosmicRayShieldGeom_CRSAbsorberLayer_hh */
