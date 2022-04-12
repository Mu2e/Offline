#ifndef PTMGeom_PTMStand_hh
#define PTMGeom_PTMStand_hh

#include <memory>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/GeomPrimitives/inc/ExtrudedSolid.hh"
#include "Offline/GeomPrimitives/inc/Box.hh"

// ProductionTarget Monitor Stand (PTM) 
//
// Author: Helenka Casler
//

// TODO: ADD THE FOLLOWING:
// - aluminum base that splays out from the center column
// - center column that keeps it at the right height
// - support blocks on the center column
// - top wedge that keeps the detectors at the right angle

namespace mu2e {

  class PTMStand {

  public:
    PTMStand(std::shared_ptr<ExtrudedSolid> topWedge,
          std::shared_ptr<Box> wedgeCutout,
          CLHEP::Hep3Vector wedgeCutoutRelPosition,
          CLHEP::HepRotation wedgeCutoutRelRotation,
          std::shared_ptr<Box> columnExtrusion,
          std::vector<CLHEP::Hep3Vector> columnOriginsInMu2e,
          std::string columnMateriaName,
          std::string wedgeMaterialName);
    PTMStand() {}

    const GenericTrap* topWedge() const { return _topWedge.get(); }
    const Box* wedgeCutout() const { return _wedgeCutout.get(); }
    CLHEP::Hep3Vector const & wedgeCutoutRelPosition() const { return _wedgeCutoutRelPosition; }
    CLHEP::HepRotation const & wedgeCutoutRelRotation() const { return _wedgeCutoutRelRotation; }
    const Box* columnExtrusion() const {return _columnExtrusion.get(); }
    const std::vector<CLHEP::Hep3Vector> columnOriginsInMu2e()   const { return _columnOriginsInMu2e; }
    const std::string columnMateriaName() const { return _columnMaterialName; }
    const std::string wedgeMaterialName() const { return _wedgeMaterialName; }




  private:
    std::shared_ptr<ExtrudedSolid> _topWedge; // contains its own origin in Mu2e and material name
    std::shared_ptr<Box> _wedgeCutout;
    CLHEP::Hep3Vector _wedgeCutoutRelPosition;
    CLHEP::HepRotation _wedgeCutoutRelRotation;
    std::shared_ptr<Box> _columnExtrusion;
    std::vector<CLHEP::Hep3Vector> _columnOriginsInMu2e;
    std::string _columnMaterialName;
    std::string _wedgeMaterialName;



  }; // class PTMStand

} // namespace mu2e

#endif