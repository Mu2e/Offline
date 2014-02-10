//
// Representation of one Absorber Layer in CosmicRayShield
//
//
// $Id: CRSAbsorberLayer.cc,v 1.1 2014/02/10 14:23:03 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2014/02/10 14:23:03 $
//
// Original author KLG based on Rob Kutschke's Layer
//

#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSAbsorberLayer.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

#ifndef __CINT__

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e 
{

  CRSAbsorberLayer::CRSAbsorberLayer():
    _id(CRSScintillatorLayerId())
  {}

  CRSAbsorberLayer::CRSAbsorberLayer(CRSScintillatorLayerId const& id):
    _id(id)
  {}

  string CRSAbsorberLayer::name( string const& base ) const
  {
    ostringstream os;
    os << base
       << _id.getShieldNumber() << "_"
       << _id.getModuleNumber() << "_"
       << _id.getLayerNumber();
    return os.str();
  }

} // namespace mu2e
#endif

