//
// Representation of one aluminum sheet in CosmicRayShield
//
//
// $Id: CRSAbsorberLayer.cc,v 1.1 2014/02/10 14:23:03 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/10 14:23:03 $
//

#include "CosmicRayShieldGeom/inc/CRSAluminumSheet.hh"

namespace mu2e 
{

  CRSAluminumSheet::CRSAluminumSheet(const CLHEP::Hep3Vector &position, const std::vector<double> &halfLength) : 
  _position(position),
  _halfLengths(halfLength)
  {}

} // namespace mu2e
