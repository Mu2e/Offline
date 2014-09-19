// Create heat and radiation shield inside the PS.
//
// Original author Andrei Gaponenko

#include <iostream>

#include "Mu2eG4/inc/constructPSShield.hh"

#include "cetlib/exception.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4Color.hh"

#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "G4Helper/inc/VolumeInfo.hh"

#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/finishNesting.hh"

#include "ProductionSolenoidGeom/inc/PSShield.hh"
#include "ProductionSolenoidGeom/inc/ProductionSolenoid.hh"

#include "GeometryService/inc/G4GeometryOptions.hh"
#include "GeometryService/inc/GeometryService.hh"
#include "GeometryService/inc/GeomHandle.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {
  void constructPSShield(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<PSShield> hrs;
    GeomHandle<ProductionSolenoid> ps;

    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "PSShield", "PSShield" );

    // Do not want to put these lookups in a loop, since it would
    // increase the processing time
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible( "PSShield" );
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck     ( "PSShield" );
    const bool placePV             = geomOptions->placePV            ( "PSShield" );
    const bool isVisible           = geomOptions->isVisible          ( "PSShield" );
    const bool isSolid             = geomOptions->isSolid            ( "PSShield" );

    // place the polycone shells
    for(unsigned ishell=0; ishell < hrs->shells().size(); ++ishell) {
      std::ostringstream osnum;
      osnum << ishell + 1;

      const Polycone& shell = hrs->shells()[ishell];
      G4VSolid *psssolid = reg.add(new G4Polycone("PSShieldPoly"+osnum.str(),
                                                  0, 2*M_PI,
                                                  shell.numZPlanes(),
                                                  &shell.zPlanes()[0],
                                                  &shell.rInner()[0],
                                                  &shell.rOuter()[0]
                                                  )
                                   );

      VolumeInfo pss("PSShieldShell"+osnum.str(),
                     shell.originInMu2e() - parent.centerInMu2e(),
                     parent.centerInWorld);

      pss.solid = psssolid;
      pss.solid->SetName(pss.name);

      finishNesting(pss,
                    findMaterialOrThrow(shell.materialName()),
                    0,
                    pss.centerInParent,
                    parent.logical,
                    0,
		    isVisible,
                    //G4Colour(0xFF/double(0xFF), 0x99/double(0xFF), 0),
                    G4Colour(config.getHep3Vector("PSShield.color"+osnum.str())),
		    isSolid,
                    forceAuxEdgeVisible,
                    placePV,
                    doSurfaceCheck
                    );
    }

//FIXME: *** the following piece of code needs to become
//FIXME: *** a double loop over shells and grooves
//FIXME:
//FIXME:    // and cut out the grooves
//FIXME:    for(unsigned i=0; i<hrs->nGrooves(); ++i) {
//FIXME:      std::ostringstream gname;
//FIXME:      gname<<"PSShieldGroove"<<i;
//FIXME:      G4Tubs *groove = reg.add(new G4Tubs(gname.str(),
//FIXME:                                          0., hrs->grooves()[i].r(),
//FIXME:                                          hrs->grooves()[i].halfLength(),
//FIXME:                                          0., 2*M_PI
//FIXME:                                          ));
//FIXME:
//FIXME:      std::ostringstream tmpname;
//FIXME:      tmpname<<"PSShieldTemp"<<i;
//FIXME:      psssolid = new G4SubtractionSolid(tmpname.str(),
//FIXME:                                        psssolid,
//FIXME:                                        groove,
//FIXME:                                        hrs->grooves()[i].placement()
//FIXME:                                        );
//FIXME:
//FIXME:    }
    if(hrs->nGrooves() > 0) {
      throw cet::exception("GEOM")<<"in constructPSShield(): only nGrooves==0 is implemented at the moment.  Can't build the requested geometry with nGrooves = "<<hrs->nGrooves()<<"\n";
    }

  } // end Mu2eWorld::constructPSShield
}
