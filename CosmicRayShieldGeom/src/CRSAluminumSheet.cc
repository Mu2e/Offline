//
// Representation of one aluminum sheet in CosmicRayShield
//
//
//

#include "Offline/CosmicRayShieldGeom/inc/CRSAluminumSheet.hh"

namespace mu2e
{

  CRSAluminumSheet::CRSAluminumSheet(const CLHEP::Hep3Vector &position, const std::vector<double> &halfLength) :
  _position(position),
  _halfLengths(halfLength)
  {}

} // namespace mu2e
