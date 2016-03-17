#ifndef CosmicRayShieldGeom_CRSSupportStructure_hh
#define CosmicRayShieldGeom_CRSSupportStructure_hh
//
// Representation of the support structure in  CosmicRayShield
//
// $Id: CRSAbsorberLayer.hh,v 1.1 2014/02/10 14:23:03 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/10 14:23:03 $
//

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e 
{

  class CRSSupportStructure
  {

    friend class CRSScintillatorModule;
    friend class CRSScintillatorShield;
    friend class CosmicRayShieldMaker;

    public:

    CRSSupportStructure() {}

    CRSSupportStructure(const std::string &name, const CLHEP::Hep3Vector &position, const std::vector<double> &halfLengths, const std::string &materialName);

    // Accept the compiler generated destructor, copy constructor and assignment operators

    const std::string &getName() const                {return _name;}
    const CLHEP::Hep3Vector &getPosition() const      {return _position;}
    const std::vector<double> &getHalfLengths() const {return _halfLengths;}
    const std::string &getMaterialName() const        {return _materialName;}

    private:

    std::string              _name;
    CLHEP::Hep3Vector        _position;
    std::vector<double>      _halfLengths;
    std::string              _materialName;
  };
}

#endif /* CosmicRayShieldGeom_CRSSupportStructure_hh */
