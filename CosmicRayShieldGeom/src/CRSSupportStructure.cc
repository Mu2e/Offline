//
// Representation of one FEB in CosmicRayShield
//
//
//

#include "CosmicRayShieldGeom/inc/CRSSupportStructure.hh"

namespace mu2e 
{

  CRSSupportStructure::CRSSupportStructure(const std::string &name, const CLHEP::Hep3Vector &position, const std::vector<double> &halfLengths, const std::string &materialName) :
  _name(name),
  _position(position),
  _halfLengths(halfLengths),
  _materialName(materialName)
  {}

} // namespace mu2e
