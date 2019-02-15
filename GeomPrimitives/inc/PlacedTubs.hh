#ifndef GeometryPrimitives_PlacedTubs_hh
#define GeometryPrimitives_PlacedTubs_hh
//
//  Properties needed to create a Geant4 tube section object:
//   - Shape parameters of the TUBS
//   - Placement information: position and rotation
//   - Material name
//
//  $Id: PlacedTubs.hh,v 1.1 2013/01/07 03:55:10 kutschke Exp $
//  $Author: kutschke $
//  $Date: 2013/01/07 03:55:10 $
//
//  Original author Rob Kutschke
//

#include <string>
#include <ostream>

#include "GeomPrimitives/inc/TubsParams.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

namespace mu2e {

  class PlacedTubs{

  public:

    PlacedTubs();

    PlacedTubs( std::string const&        name,
                TubsParams const&         params,
                CLHEP::Hep3Vector const&  position,
                CLHEP::HepRotation const& rotation,
                std::string const&        materialName);

    PlacedTubs(  std::string const&       name,
                 TubsParams const&        params,
                 CLHEP::Hep3Vector const& position,
                 std::string const&       materialName);

    // Accept compiler supplied destructor, copy c'tor and assignment operator.


    std::string        const& name()         const { return _name; }
    TubsParams         const& tubsParams()   const { return _params;  }
    CLHEP::Hep3Vector  const& position()     const { return _position; }
    CLHEP::HepRotation const& rotation()     const { return _rotation; }
    std::string        const& materialName() const { return _materialName; }


    // Forward the accessors of TubsParams.
    double innerRadius()  const { return _params.innerRadius(); }
    double outerRadius()  const { return _params.outerRadius(); }
    double zHalfLength()  const { return _params.zHalfLength(); }
    double phi0()         const { return _params.phi0();        }
    double phiMax()       const { return _params.phiMax();      }
    double const * data() const { return _params.data();        }

  private:

    std::string        _name;          // Name used for G4 volume names
    TubsParams         _params;
    CLHEP::Hep3Vector  _position;      // Relative to mother volume.
    CLHEP::HepRotation _rotation;      // Relative to mother volume.
    std::string        _materialName;

  };

  inline std::ostream& operator<<(std::ostream& ost,
                                  const PlacedTubs& pt ){
    ost << "name:        " << pt.name()         << "\n"
        << "Tubs Params: " << pt.tubsParams()   << "\n"
        << "Position:    " << pt.position()     << "\n"
        << "Rotation:    " << pt.rotation()     << "\n"
        << "Material:    " << pt.materialName();
    return ost;
  }

}

#endif /* GeomPrimitives_PlacedTubs_hh */
