//
// Original author Andrei Gaponenko

#include <sstream>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "GeometryService/inc/PSVacuumMaker.hh"
#include "ProductionSolenoidGeom/inc/PSVacuum.hh"

#include "ProductionSolenoidGeom/inc/ProductionSolenoid.hh"
#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

namespace mu2e {

  std::unique_ptr<PSVacuum>  PSVacuumMaker::make(const SimpleConfig& c,
                                               const ProductionSolenoid& ps,
                                               const PSEnclosure& pse,
                                               double zmax
                                               )

  {
    double pseHL = pse.shell().halfLength();
    if (pseHL < 1.0e-06) pseHL = pse.shellCone().halfLength();
    const double zmin = ps.psEndRefPoint().z() - 2*pseHL + pse.getExtraOffset();

    const CLHEP::Hep3Vector originInMu2e(ps.psEndRefPoint().x(),
                                         ps.psEndRefPoint().y(),
                                         (zmin+zmax)/2);

    std::unique_ptr<PSVacuum> psv(new PSVacuum(
                                               Tube(
                                                    c.getString("PS.vacuumMaterialName"),
                                                    originInMu2e,
                                                    0.,
                                                    ps.getVacVesselInnerParamsPtr()->innerRadius(),
                                                    (zmax - zmin)/2
                                                    )
                                               )
                                  );


    psv->_vacuumPressure   = c.getDouble("PS.vacuumPressure");
    psv->_vacuumG4Material = c.getString("PS.vacuumG4Material");

    return psv;

  }

} // namespace mu2e
