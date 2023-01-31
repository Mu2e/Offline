#include "Offline/PTMGeom/inc/PTMStand.hh"

// ProductionTarget Monitor Stand (PTM) Object
//
// Author: Helenka Casler
// add a comment

namespace mu2e{
  PTMStand::PTMStand(std::shared_ptr<ExtrudedSolid> topWedge,
          double wedgeCutoutSemiMajor,
          double wedgeCutoutSemiMinor,
          CLHEP::Hep3Vector wedgeCutoutRelPosition,
          CLHEP::HepRotation wedgeCutoutRelRotation,
          std::shared_ptr<Box> columnExtrusion,
          std::vector<CLHEP::Hep3Vector> columnOriginsInMu2e,
          std::vector<CLHEP::HepRotation> columnRotations,
          std::string columnMaterialName,
          std::string wedgeMaterialName)
          : _topWedge(topWedge),
            _wedgeCutoutSemiMajor(wedgeCutoutSemiMajor),
            _wedgeCutoutSemiMinor(wedgeCutoutSemiMinor),
            _wedgeCutoutRelPosition(wedgeCutoutRelPosition),
            _wedgeCutoutRelRotation(wedgeCutoutRelRotation),
            _columnExtrusion(columnExtrusion),
            _columnOriginsInMu2e(columnOriginsInMu2e),
            _columnRotations(columnRotations),
            _columnMaterialName(columnMaterialName),
            _wedgeMaterialName(wedgeMaterialName) {

  }

}
