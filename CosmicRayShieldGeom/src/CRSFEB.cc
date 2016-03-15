//
// Representation of one FEB in CosmicRayShield
//
//
// $Id: CRSAbsorberLayer.cc,v 1.1 2014/02/10 14:23:03 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/10 14:23:03 $
//

#include "CosmicRayShieldGeom/inc/CRSFEB.hh"

namespace mu2e 
{

  CRSFEB::CRSFEB(const CLHEP::Hep3Vector &position, const std::vector<double> &halfLength) : 
  _position(position),
  _halfLengths(halfLength)
  {}

} // namespace mu2e
