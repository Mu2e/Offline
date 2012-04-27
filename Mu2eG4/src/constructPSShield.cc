// Create heat and radiation shield inside the PS.
//
// Original author Andrei Gaponenko

#include "Mu2eG4/inc/constructPSShield.hh"

#include "CLHEP/Vector/ThreeVector.h"

#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Color.hh"

#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "G4Helper/inc/VolumeInfo.hh"

#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/finishNesting.hh"

#include "ProductionSolenoidGeom/inc/PSShield.hh"
#include "ProductionSolenoidGeom/inc/ProductionSolenoid.hh"
#include "GeometryService/inc/GeomHandle.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {
  void constructPSShield(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<PSShield> hrs;
    GeomHandle<ProductionSolenoid> ps;

    MaterialFinder materialFinder(config);
    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    const bool placePV             = true;

    // start with the polycone
    G4VSolid *psssolid = reg.add(new G4Polycone("PSShieldPoly",
                                                0, 2*M_PI,
                                                hrs->bulk().numZPlanes(),
                                                &hrs->bulk().zPlanes()[0],
                                                &hrs->bulk().rInner()[0],
                                                &hrs->bulk().rOuter()[0]
                                                )
                                 );

    // and cut out the grooves
    for(unsigned i=0; i<hrs->nGrooves(); ++i) {
      std::ostringstream gname;
      gname<<"PSShieldGroove"<<i;
      G4Tubs *groove = reg.add(new G4Tubs(gname.str(),
                                          0., hrs->grooves()[i].r(),
                                          hrs->grooves()[i].halfLength(),
                                          0., 2*M_PI
                                          ));

      std::ostringstream tmpname;
      tmpname<<"PSShieldTemp"<<i;
      psssolid = new G4SubtractionSolid(tmpname.str(),
                                        psssolid,
                                        groove,
                                        hrs->grooves()[i].placement()
                                        );

    }

    //----------------------------------------------------------------
    VolumeInfo pss("PSShield", hrs->originInMu2e() - parent.centerInMu2e(), parent.centerInWorld);
    pss.solid = psssolid;
    pss.solid->SetName(pss.name);

    finishNesting(pss,
                  materialFinder.get("PSShield.materialName"),
                  0,
                  pss.centerInParent,
                  parent.logical,
                  0,
                  config.getBool("PSShield.visible"),
                  G4Colour(0xFF/double(0xFF), 0x99/double(0xFF), 0),
                  config.getBool("PSShield.solid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

  } // end Mu2eWorld::constructPSShield
}
