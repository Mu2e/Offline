#include "ProductionTargetGeom/inc/ProductionTargetMaker.hh"

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

#include "ProductionTargetGeom/inc/ProductionTarget.hh"

namespace mu2e {

  ProductionTargetMaker::ProductionTargetMaker(const SimpleConfig& c, double solenoidOffset)
    : m_det(new ProductionTarget(c.getDouble("targetPS_rOut"),
                                 c.getDouble("targetPS_halfLength"),
                                 c.getDouble("targetPS_rotX") * CLHEP::degree,
                                 c.getDouble("targetPS_rotY") * CLHEP::degree,
                                 CLHEP::Hep3Vector(solenoidOffset,
                                                   0,
                                                   c.getDouble("productionTarget.zNominal")
                                                   )
                                 + c.getHep3Vector("productionTarget.offset", CLHEP::Hep3Vector(0,0,0))
                                 )
            )
  {}

}
