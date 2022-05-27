//
// Representation of one Absorber Layer in CosmicRayShield
//
//
//

#include "Offline/CosmicRayShieldGeom/inc/CRSAbsorberLayer.hh"

namespace mu2e
{

  CRSAbsorberLayer::CRSAbsorberLayer(const CLHEP::Hep3Vector &position, const std::vector<double> &halfLength) :
  _position(position),
  _halfLengths(halfLength)
  {}

} // namespace mu2e
