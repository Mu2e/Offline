//
// Representation of one Scintillator Layer in CosmicRayShield
//
//
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

  CRSScintillatorLayer::CRSScintillatorLayer()
  {
    _halfLengths.resize(3);
    _localToWorld.resize(3);
  }

  CRSScintillatorLayer::CRSScintillatorLayer(CRSScintillatorLayerId const& id) :
    _id(id)
  {
    _halfLengths.resize(3);
    _localToWorld.resize(3);
  }

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

