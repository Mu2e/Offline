#ifndef CosmicRayShieldGeom_CRSAluminumSheet_hh
#define CosmicRayShieldGeom_CRSAluminumSheet_hh
//
// Representation of one aluminum sheet in  CosmicRayShield
//
// $Id: CRSAbsorberLayer.hh,v 1.1 2014/02/10 14:23:03 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/10 14:23:03 $
//

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e 
{

  class CRSAluminumSheet
  {

    friend class CRSScintillatorModule;
    friend class CRSScintillatorShield;
    friend class CosmicRayShieldMaker;

    public:

    CRSAluminumSheet() {}

    CRSAluminumSheet(const CLHEP::Hep3Vector &position, const std::vector<double> &halfLength);
    // Accept the compiler generated destructor, copy constructor and assignment operators

    const CLHEP::Hep3Vector &getPosition() const {return _position;}
    const std::vector<double> &getHalfLengths() const {return _halfLengths;}

    private:

    CLHEP::Hep3Vector _position;
    std::vector<double> _halfLengths;
  };
}

#endif /* CosmicRayShieldGeom_CRSAluminumSheet_hh */
