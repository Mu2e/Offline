//
// Representation of one FEB in CosmicRayShield
//
//
// $Id: CRSAbsorberLayer.cc,v 1.1 2014/02/10 14:23:03 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/10 14:23:03 $
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
