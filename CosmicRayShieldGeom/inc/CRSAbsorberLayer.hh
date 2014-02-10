#ifndef CosmicRayShieldGeom_CRSAbsorberLayer_hh
#define CosmicRayShieldGeom_CRSAbsorberLayer_hh
//
// Representation of one Absorber Layer in  CosmicRayShield
//
// $Id: CRSAbsorberLayer.hh,v 1.1 2014/02/10 14:23:03 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/10 14:23:03 $
//
// Original author KLG; somewhat based on  Rob Kutschke's Layer
//

#include <vector>
#include <deque>

#include "CosmicRayShieldGeom/inc/CRSScintillatorLayerId.hh"

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e 
{

  class CRSAbsorberLayer
  {

    friend class CRSScintillatorModule;
    friend class CRSScintillatorShield;
    friend class CosmicRayShieldMaker;

    public:

    CRSAbsorberLayer();

    CRSAbsorberLayer(CRSScintillatorLayerId const & id);

    // Accept the compiler generated destructor, copy constructor and assignment operators

    CRSScintillatorLayerId const & id() const { return _id;}

    const CLHEP::Hep3Vector &getPosition() const {return _position;}
    const std::vector<double> &getHalfLengths() const {return _halfLengths;}

    // Formatted string embedding the id of the layer.
    std::string name( std::string const & base ) const;

    private:

    CRSScintillatorLayerId _id;

    CLHEP::Hep3Vector _position;
    std::vector<double> _halfLengths;
  };
}

#endif /* CosmicRayShieldGeom_CRSAbsorberLayer_hh */
