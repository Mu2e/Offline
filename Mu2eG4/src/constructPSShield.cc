// Create heat and radiation shield inside the PS.
//
// Original author Andrei Gaponenko

#include <iostream>

#include "Mu2eG4/inc/constructPSShield.hh"

#include "cetlib_except/exception.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Geant4/G4Tubs.hh"
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4UnionSolid.hh"
#include "Geant4/G4Color.hh"

#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/nestPolycone.hh"
#include "Mu2eG4/inc/nestTubs.hh"

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

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "PSShield", "PSShield" );

    // Do not want to put these lookups in a loop, since it would
    // increase the processing time
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible( "PSShield" );
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck     ( "PSShield" );
    const bool isVisible           = geomOptions->isVisible          ( "PSShield" );
    const bool isSolid             = geomOptions->isSolid            ( "PSShield" );
    const bool placePV             = geomOptions->placePV            ( "PSShield" );


    const int version              = hrs->version();

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    // -----------------------------
    // Put in beam pipe inlet.  David Norvil Brown, Louisville, March 2015

    Tube const & pssInletParams = hrs->beamInlet();
    //    G4Material* beampipeMaterial = findMaterialOrThrow(pssInletParams.materialName());
    CLHEP::Hep3Vector place = hrs->getBeamInletCenter();
    CLHEP::HepRotation * turn = reg.add(CLHEP::HepRotation(CLHEP::HepRotation::IDENTITY));
    turn->rotateY(hrs->getBeamAngleY());
    turn->rotateX(hrs->getBeamAngleX());

    // Keeps overlapping with PSVacuum
    //    VolumeInfo pssBeamInletInfo =
    //      nestTubs("PSSBeamInlet", pssInletParams.getTubsParams(),
    //                 beampipeMaterial, turn,  place-parent.centerInMu2e(),
    //                 parent, 0, G4Colour::Blue(), "PSShield");

    // Now make the volume to subtract from HRS.
    G4Tubs* beamPassTub = new G4Tubs( "beampipePassthrough",
                                      0.0,
                                      pssInletParams.getTubsParams().data()[1]+2*CLHEP::mm,
                                      pssInletParams.getTubsParams().data()[2]+600.0*CLHEP::mm,
                                      0.0*CLHEP::degree,
                                      360.0*CLHEP::degree );


    // ------------------------------------------------------------
    // place the polycone shells for the HRS proper
    for(unsigned ishell=0; ishell < hrs->shells().size(); ++ishell) {
      std::ostringstream osnum;
      osnum << ishell + 1;

      const Polycone& shell = hrs->shells()[ishell];
      G4VSolid *psssolid = new G4Polycone("PSShieldPoly"+osnum.str(),
                                          0, 2*M_PI,
                                          shell.numZPlanes(),
                                          &shell.zPlanes()[0],
                                          &shell.rInner()[0],
                                          &shell.rOuter()[0]
                                          );

      VolumeInfo pss("PSShieldShell"+osnum.str(),
                     shell.originInMu2e() - parent.centerInMu2e(),
                     parent.centerInWorld);


      //      std::cout << "inside " << __func__ << "PSSShieldShell"+osnum.str() << " " << shell.originInMu2e() << " " << parent.centerInMu2e() << parent.centerInWorld << std::endl;
      G4SubtractionSolid* aSolid =
        new G4SubtractionSolid ( pss.name,
                                 psssolid,
                                 beamPassTub,
                                 turn,
                                 place-shell.originInMu2e());

      // pss.solid = psssolid;
      pss.solid = aSolid;
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
      //    std::cout << "and after finish nesting = " << pss.centerInParent << std::endl;
    }
    // Insert the endRings if this is version 2 or higher
    // Treat the extra water and stainless sheath as end rings.
    if ( version > 1 ) {
      for(unsigned iER=0; iER < hrs->endRings().size(); ++iER) {
        std::ostringstream osnum;
        osnum << iER + 1;

        const Polycone& eRing = hrs->endRings()[iER];
        G4VSolid *psersolid = new G4Polycone("PSShieldEndRing"+osnum.str(),
                                            0, 2*M_PI,
                                            eRing.numZPlanes(),
                                            &eRing.zPlanes()[0],
                                            &eRing.rInner()[0],
                                            &eRing.rOuter()[0]
                                            );

        VolumeInfo pser("PSShieldEndRing"+osnum.str(),
                       eRing.originInMu2e() - parent.centerInMu2e(),
                       parent.centerInWorld);
        //      std::cout << "endring = " << eRing.originInMu2e() << " " << parent.centerInMu2e() << std::endl;
        G4VSolid *aSolid = 0;
        if ( 1 == iER ) {
          // downstream - beam inlet tube penetrates
          aSolid = new G4SubtractionSolid ( pser.name,
                                            psersolid,
                                            beamPassTub,
                                            turn,
                                            place-eRing.originInMu2e());
        } else {
          // upstream - beam inlet tube does not penetrate
          // Same for extra water and stainless sheath.
          aSolid = psersolid;
        }

        pser.solid = aSolid;
        pser.solid->SetName(pser.name);

        finishNesting(pser,
                      findMaterialOrThrow(eRing.materialName()),
                      0,
                      pser.centerInParent,
                      parent.logical,
                      0,
                      isVisible,
                      G4Colour(config.getHep3Vector("PSShield.endRing.color")),
                      isSolid,
                      forceAuxEdgeVisible,
                      placePV,
                      doSurfaceCheck
                      );
      }

    } // end if ( version > 1)...  (putting in endRings)


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
