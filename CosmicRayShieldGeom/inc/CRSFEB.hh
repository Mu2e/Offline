#ifndef CosmicRayShieldGeom_CRSFEB_hh
#define CosmicRayShieldGeom_CRSFEB_hh
//
// Representation of one FEB in  CosmicRayShield
//
//

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e 
{

  class CRSFEB
  {

    friend class CRSScintillatorModule;
    friend class CRSScintillatorShield;
    friend class CosmicRayShieldMaker;

    public:

    CRSFEB() {}

    CRSFEB(const CLHEP::Hep3Vector &position, const std::vector<double> &halfLength);

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const CLHEP::Hep3Vector &getPosition() const {return _position;}
    const std::vector<double> &getHalfLengths() const {return _halfLengths;}

    private:

    CLHEP::Hep3Vector _position;
    std::vector<double> _halfLengths;
  };
}

#endif /* CosmicRayShieldGeom_CRSFEB_hh */
