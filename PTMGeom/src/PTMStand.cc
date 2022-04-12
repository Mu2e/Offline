#include "Offline/PTMGeom/inc/PTMStand.hh"

// ProductionTarget Monitor Stand (PTM) Object
//
// Author: Helenka Casler
//

namespace mu2e{
  PTMStand::PTMStand(std::shared_ptr<ExtrudedSolid> topWedge,
          std::shared_ptr<Box> wedgeCutout,
          CLHEP::Hep3Vector wedgeCutoutRelPosition,
          CLHEP::HepRotation wedgeCutoutRelRotation,
          std::shared_ptr<Box> columnExtrusion,
          std::vector<CLHEP::Hep3Vector> columnOriginsInMu2e, 
          std::vector<CLHEP::HepRotation> columnRotationsInMu2e,
          std::string columnMateriaName,
          std::string wedgeMaterialName)
          : _topWedge(topWedge),
            _wedgeCutout(wedgeCutout),
            _wedgeCutoutRelPosition(wedgeCutoutRelPosition),
            _wedgeCutoutRelRotation(wedgeCutoutRelRotation),
            _columnExtrusion(columnExtrusion),
            _columnOriginsInMu2e(columnOriginsInMu2e), 
            _columnRotationsInMu2e(columnRotationsInMu2e),
            _columnMateriaName(columnMateriaName),
            _wedgeMaterialName(wedgeMaterialName) {

  }

}
