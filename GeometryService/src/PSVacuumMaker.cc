//
// Original author Andrei Gaponenko

#include <cmath>
#include <sstream>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "Offline/GeometryService/inc/PSVacuumMaker.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSVacuum.hh"

#include "Offline/ProductionSolenoidGeom/inc/ProductionSolenoid.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSEnclosure.hh"

#include "Offline/BeamlineGeom/inc/TransportSolenoid.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"

#include "cetlib_except/exception.h"

namespace mu2e {

  std::unique_ptr<PSVacuum>  PSVacuumMaker::make(const SimpleConfig& c,
                                               const ProductionSolenoid& ps,
                                               const PSEnclosure& pse,
                                               double zmax,
                                               const TransportSolenoid& ts
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

    // Conical hard connect between the PS bore tube and the TSu ring plate
    // (Mu2e-doc-8545). It starts at the TS end of the PS vacuum vessel inner
    // tube, with that tube's inner radius, and ends on the upstream face of the
    // first TS ring, which must be the ring around TS1 on the PS axis.
    if (c.getBool("PS.TSInterconnect.build", false)) {
      const Tube& inner = *ps.getVacVesselInnerParamsPtr();
      const double zPS = inner.originInMu2e().z() + inner.halfLength();

      if (ts.xRing().empty() || std::abs(ts.thetaRing()[0]) > 1.e-9
          || std::abs(ts.xRing()[0] - inner.originInMu2e().x()) > 1.e-6
          || std::abs(ts.yRing()[0] - inner.originInMu2e().y()) > 1.e-6) {
        throw cet::exception("GEOM")
          << "PSVacuumMaker: PS.TSInterconnect.build needs TS ring 0 on the PS axis with ts.thetaRing = 0\n";
      }
      const double zTS = ts.zRing()[0] - ts.lengthRing()/2. - ts.thickRingSide();
      const double length = zTS - zPS;
      if (length <= 0.) {
        throw cet::exception("GEOM")
          << "PSVacuumMaker: TS ring 0 face at z = " << zTS
          << " is not downstream of the PS inner vessel end at z = " << zPS << "\n";
      }

      const double rIn1 = inner.innerRadius();
      const double rIn2 = c.getDouble("PS.TSInterconnect.rIn2");
      // Radial thickness that gives the requested wall thickness normal to the cone surface
      const double cosHalfAngle = length/std::hypot(length, rIn2 - rIn1);
      const double radialThickness = c.getDouble("PS.TSInterconnect.thickness")/cosHalfAngle;
      // The cone end must face the ring plate, at least in part
      if (rIn2 >= ts.rOutRingSide() || rIn2 + radialThickness <= ts.rInRingSide()) {
        throw cet::exception("GEOM")
          << "PSVacuumMaker: the PS-TS cone end, r = " << rIn2 << " to " << rIn2 + radialThickness
          << ", does not face the TS ring plate, r = " << ts.rInRingSide()
          << " to " << ts.rOutRingSide() << "\n";
      }

      psv->tsInterconnect_ = Cone(rIn1, rIn1 + radialThickness,
                                  rIn2, rIn2 + radialThickness,
                                  length/2., 0., 2.*M_PI,
                                  CLHEP::Hep3Vector(inner.originInMu2e().x(),
                                                    inner.originInMu2e().y(),
                                                    (zPS + zTS)/2.),
                                  CLHEP::HepRotation(),
                                  c.getString("PS.TSInterconnect.materialName"));
      psv->hasTSInterconnect_ = true;
    }

    return psv;

  }

} // namespace mu2e
