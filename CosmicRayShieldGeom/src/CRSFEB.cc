//
// Representation of one FEB in CosmicRayShield
//
//
//

#include "Offline/CosmicRayShieldGeom/inc/CRSFEB.hh"

namespace mu2e
{

  CRSFEB::CRSFEB(const CLHEP::Hep3Vector &position, const std::vector<double> &halfLength) :
  _position(position),
  _halfLengths(halfLength)
  {}

} // namespace mu2e
