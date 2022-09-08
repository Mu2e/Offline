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

namespace mu2e {

  class PTMStand {

  public:
    PTMStand(std::shared_ptr<ExtrudedSolid> topWedge,
          double wedgeCutoutSemiMajor,
          double wedgeCutoutSemiMinor,
          CLHEP::Hep3Vector wedgeCutoutRelPosition,
          CLHEP::HepRotation wedgeCutoutRelRotation,
          std::shared_ptr<Box> columnExtrusion,
          std::vector<CLHEP::Hep3Vector> columnOriginsInMu2e,
          std::vector<CLHEP::HepRotation> columnRotations,
          std::string columnMaterialName,
          std::string wedgeMaterialName);
    PTMStand() {}

    const ExtrudedSolid* topWedge() const { return _topWedge.get(); }
    double wedgeCutoutSemiMajor() const { return _wedgeCutoutSemiMajor; }
    double wedgeCutoutSemiMinor() const { return _wedgeCutoutSemiMinor; }
    CLHEP::Hep3Vector const & wedgeCutoutRelPosition() const { return _wedgeCutoutRelPosition; }
    CLHEP::HepRotation const & wedgeCutoutRelRotation() const { return _wedgeCutoutRelRotation; }
    const Box* columnExtrusion() const {return _columnExtrusion.get(); }
    const std::vector<CLHEP::Hep3Vector> columnOriginsInMu2e()   const { return _columnOriginsInMu2e; }
    const std::vector<CLHEP::HepRotation> columnRotations()   const { return _columnRotations; }
    const std::string columnMaterialName() const { return _columnMaterialName; }
    const std::string wedgeMaterialName() const { return _wedgeMaterialName; }




  private:
    std::shared_ptr<ExtrudedSolid> _topWedge; // contains its own origin in Mu2e and material name
    double _wedgeCutoutSemiMajor;
    double _wedgeCutoutSemiMinor;
    CLHEP::Hep3Vector _wedgeCutoutRelPosition;
    CLHEP::HepRotation _wedgeCutoutRelRotation;
    std::shared_ptr<Box> _columnExtrusion;
    std::vector<CLHEP::Hep3Vector> _columnOriginsInMu2e;
    std::vector<CLHEP::HepRotation> _columnRotations;
    std::string _columnMaterialName;
    std::string _wedgeMaterialName;



  }; // class PTMStand

} // namespace mu2e

#endif
