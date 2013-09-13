//
// Representation of one Scintillator Layer in CosmicRayShield
//
//
// $Id: CRSScintillatorLayer.cc,v 1.4 2013/09/13 06:42:44 ehrlich Exp $
// $Author: ehrlich $
// $Date: 2013/09/13 06:42:44 $
//
// Original author KLG based on Rob Kutschke's Layer
//

#include <sstream>

#include "CosmicRayShieldGeom/inc/CRSScintillatorLayer.hh"
#include "CosmicRayShieldGeom/inc/CosmicRayShield.hh"

#ifndef __CINT__

using namespace std;

using CLHEP::Hep3Vector;

namespace mu2e 
{

  CRSScintillatorLayer::CRSScintillatorLayer():
    _id(CRSScintillatorLayerId())
  {}

  CRSScintillatorLayer::CRSScintillatorLayer(CRSScintillatorLayerId const& id):
    _id(id)
  {}

  string CRSScintillatorLayer::name( string const& base ) const
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

