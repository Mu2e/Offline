//
//  Properties needed to create a Geant4 tube section object:
//   - Shape parameters of the TUBS
//   - Placement information: position and rotation
//   - Material name
//
//
//  Original author Rob Kutschke
//

#include "GeomPrimitives/inc/PlacedTubs.hh"

namespace mu2e {

  PlacedTubs::PlacedTubs():
    _name(),
    _params(0.,0.,0.),
    _position(),
    _rotation(),
    _materialName(){
  }

  PlacedTubs::PlacedTubs( std::string const&        name,
                          TubsParams const&         params,
                          CLHEP::Hep3Vector const&  position,
                          CLHEP::HepRotation const& rotation,
                          std::string const&        materialName):
    _name(name),
    _params(params),
    _position(position),
    _rotation(rotation),
    _materialName(materialName){
  }

  PlacedTubs::PlacedTubs( std::string const&       name,
                          TubsParams const&        params,
                          CLHEP::Hep3Vector const& position,
                          std::string const&       materialName):
    _name(name),
    _params(params),
    _position(position),
    _rotation(),
    _materialName(materialName){
  }


}
