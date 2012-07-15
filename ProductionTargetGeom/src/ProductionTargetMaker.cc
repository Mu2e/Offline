#include "ProductionTargetGeom/inc/ProductionTargetMaker.hh"

#include "cetlib/exception.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "ConfigTools/inc/SimpleConfig.hh"

#include "ProductionTargetGeom/inc/ProductionTarget.hh"

namespace mu2e {

  std::auto_ptr<ProductionTarget> ProductionTargetMaker::make(const SimpleConfig& c, double solenoidOffset) {
    return std::auto_ptr<ProductionTarget>
      (new ProductionTarget(c.getDouble("targetPS_rOut"),
                            c.getDouble("targetPS_halfLength"),
                            c.getDouble("targetPS_rotX") * CLHEP::degree,
                            c.getDouble("targetPS_rotY") * CLHEP::degree,
                            CLHEP::Hep3Vector(solenoidOffset,
                                              0,
                                              c.getDouble("productionTarget.zNominal")
                                              )
                            + c.getHep3Vector("productionTarget.offset", CLHEP::Hep3Vector(0,0,0))
                            )
       );
  } // make()

} // namespace mu2e
