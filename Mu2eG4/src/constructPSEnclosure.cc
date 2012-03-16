// Andrei Gaponenko, 2012

#include "Mu2eG4/inc/constructPSEnclosure.hh"

#include "ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "GeomPrimitives/inc/Tube.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"

namespace mu2e {

  //================================================================

  void constructPSEnclosure(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<PSEnclosure> pse;

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV             = true;

    //----------------------------------------------------------------
    nestTubs("PSEnclosureShell",
             pse->shell().getTubsParams(),
             findMaterialOrThrow(pse->shell().materialName()),
             0,
             pse->shell().originInMu2e() - parent.centerInMu2e(),
             parent,
             0,
             config.getBool("PSEnclosure.visible"),
             G4Colour::Blue(),
             config.getBool("PSEnclosure.solid"),
             forceAuxEdgeVisible,
             placePV,
             doSurfaceCheck
             );

    nestTubs("PSEnclosureVacuum",
             pse->vacuum().getTubsParams(),
             findMaterialOrThrow(pse->vacuum().materialName()),
             0,
             pse->vacuum().originInMu2e() - parent.centerInMu2e(),
             parent,
             0,
             config.getBool("PSEnclosure.vacuum.visible"),
             G4Colour::Grey(),
             config.getBool("PSEnclosure.vacuum.solid"),
             forceAuxEdgeVisible,
             placePV,
             doSurfaceCheck
             );

    nestTubs("PSEnclosureEndPlate",
             pse->endPlate().getTubsParams(),
             findMaterialOrThrow(pse->endPlate().materialName()),
             0,
             pse->endPlate().originInMu2e() - parent.centerInMu2e(),
             parent,
             0,
             config.getBool("PSEnclosure.visible"),
             G4Colour::Blue(),
             config.getBool("PSEnclosure.solid"),
             forceAuxEdgeVisible,
             placePV,
             doSurfaceCheck
             );

    //----------------------------------------------------------------
  }

}
