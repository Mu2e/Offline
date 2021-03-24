// Andrei Gaponenko, 2011

#include "Mu2eG4/inc/constructProtonBeamDump.hh"
#include "Mu2eG4/inc/constructExtMonFNAL.hh"
#include "Mu2eG4Helper/inc/AntiLeakRegistry.hh"

#include <iostream>
#include <cmath>

#include "Geant4/G4Color.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4Trap.hh"
#include "Geant4/G4Orb.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4ExtrudedSolid.hh"
#include "Geant4/G4IntersectionSolid.hh"
#include "Geant4/G4TwoVector.hh"
#include "Geant4/G4SDManager.hh"

#include "Geant4/G4UniformMagField.hh"
#include "Geant4/G4Mag_UsualEqRhs.hh"
#include "Geant4/G4ExactHelixStepper.hh"
//#include "Geant4/G4NystromRK4.hh"
#include "Geant4/G4ChordFinder.hh"
#include "Geant4/G4FieldManager.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib_except/exception.h"

#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

#include "GeometryService/inc/G4GeometryOptions.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"

#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/FieldMgr.hh"

#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4Helper/inc/VolumeInfo.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  void constructProtonBeamDump(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<ProtonBeamDump> dump;
    GeomHandle<ExtMonFNALBuilding> emfb;
    GeomHandle<Mu2eEnvelope> env;
    GeomHandle<WorldG4> world;

    MaterialFinder materialFinder(config);

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "ProtonBeamDumpFront"    , "protonBeamDump.front");
    geomOptions->loadEntry( config, "ProtonBeamDumpBack"     , "protonBeamDump.back");
    geomOptions->loadEntry( config, "ProtonBeamDumpDirt"     , "protonBeamDump.dirt"       );
    geomOptions->loadEntry( config, "ProtonBeamDumpShielding", "protonBeamDump.shielding"  );
    geomOptions->loadEntry( config, "ProtonBeamDumpCore"     , "protonBeamDump.core"       );
    geomOptions->loadEntry( config, "ProtonBeamDumpMouth"    , "protonBeamDump.mouth"      );
    geomOptions->loadEntry( config, "ProtonBeamNeutronCave"  , "protonBeamDump.neutronCave");


    static CLHEP::HepRotation shieldingRot(CLHEP::HepRotation::IDENTITY);
    shieldingRot.rotateX( 90*CLHEP::degree);
    shieldingRot.rotateZ( 90*CLHEP::degree);

    VolumeInfo beamDumpFront("ProtonBeamDumpFront",
                             CLHEP::Hep3Vector(0, dump->frontShieldingCenterInMu2e()[1], 0)
                             - parent.centerInMu2e(),
                             parent.centerInWorld);

    beamDumpFront.solid = new G4ExtrudedSolid(beamDumpFront.name, dump->frontShieldingOutline(),
                                              dump->frontShieldingHalfSize()[1],
                                              G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);


    finishNesting(beamDumpFront,
                  materialFinder.get("protonBeamDump.material.shielding"),
                  &shieldingRot,
                  beamDumpFront.centerInParent,
                  parent.logical,
                  0,
                  geomOptions->isVisible( "ProtonBeamDumpFront" ),
                  G4Colour::Red() ,
                  geomOptions->isSolid( "ProtonBeamDumpFront" ),
                  geomOptions->forceAuxEdgeVisible( "ProtonBeamDumpFront" ),
                  geomOptions->placePV( "ProtonBeamDumpFront" ),
                  geomOptions->doSurfaceCheck( "ProtonBeamDumpFront" )
                  );

    VolumeInfo beamDumpBack("ProtonBeamDumpBack",
                            CLHEP::Hep3Vector(0, dump->backShieldingCenterInMu2e()[1], 0)
                            - parent.centerInMu2e(),
                            parent.centerInWorld);

    beamDumpBack.solid = new G4ExtrudedSolid(beamDumpBack.name, dump->backShieldingOutline(),
                                             dump->backShieldingHalfSize()[1],
                                             G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

    finishNesting(beamDumpBack,
                  materialFinder.get("protonBeamDump.material.shielding"),
                  &shieldingRot,
                  beamDumpBack.centerInParent,
                  parent.logical,
                  0,
                  geomOptions->isVisible( "ProtonBeamDumpBack" ),
                  G4Colour::Red() ,
                  geomOptions->isSolid( "ProtonBeamDumpBack" ),
                  geomOptions->forceAuxEdgeVisible( "ProtonBeamDumpBack" ),
                  geomOptions->placePV( "ProtonBeamDumpBack" ),
                  geomOptions->doSurfaceCheck( "ProtonBeamDumpBack" )
                  );

    static const CLHEP::HepRotation rotationInShield( (shieldingRot*dump->coreRotationInMu2e()).inverse() );
    CLHEP::Hep3Vector coreAirPositionInShield( shieldingRot * (dump->coreAirCenterInMu2e() - beamDumpFront.centerInMu2e()));
    VolumeInfo  beamDumpCoreAir = nestBox("ProtonBeamDumpCoreAir",
                                          dump->coreAirHalfSize(),
                                          materialFinder.get("protonBeamDump.material.air"),
                                          &rotationInShield,
                                          coreAirPositionInShield,
                                          beamDumpFront, 0,
                                          G4Colour::Cyan()
                                          );

    CLHEP::Hep3Vector corePositionInCoreAir( dump->coreCenterInMu2e()-dump->coreAirCenterInMu2e() );
    nestBox("ProtonBeamDumpCore",
            dump->coreHalfSize(),
            materialFinder.get("protonBeamDump.material.core"),
            0,
            corePositionInCoreAir,
            beamDumpCoreAir, 0,
            G4Colour::Blue()
            );

    CLHEP::Hep3Vector mouthPositionInShield( shieldingRot * (dump->mouthCenterInMu2e() - beamDumpFront.centerInMu2e()));
    nestBox("ProtonBeamDumpMouth",
            dump->mouthHalfSize(),
            materialFinder.get("protonBeamDump.material.air"),
            &rotationInShield,
            mouthPositionInShield,
            beamDumpFront, 0,
            G4Colour::Cyan()
            );

    CLHEP::Hep3Vector cavePositionInShield( shieldingRot * (dump->neutronCaveCenterInMu2e() - beamDumpFront.centerInMu2e()));
    nestBox("ProtonBeamNeutronCave",
            dump->neutronCaveHalfSize(),
            materialFinder.get("protonBeamDump.material.air"),
            &rotationInShield,
            cavePositionInShield,
            beamDumpFront, 0,
            G4Colour::Cyan()
            );

    CLHEP::Hep3Vector frontSteelPositionInShield( shieldingRot * (dump->frontSteelCenterInMu2e() - beamDumpFront.centerInMu2e()));
    nestBox("ProtonBeamDumpFrontSteel",
            dump->frontSteelHalfSize(),
            materialFinder.get("protonBeamDump.material.core"),
            &rotationInShield,
            frontSteelPositionInShield,
            beamDumpFront, 0,
            G4Colour::Blue()
            );

    CLHEP::Hep3Vector backSteelPositionInShield( shieldingRot * (dump->backSteelCenterInMu2e() - beamDumpBack.centerInMu2e()));
    nestBox("ProtonBeamDumpBackSteel",
            dump->backSteelHalfSize(),
            materialFinder.get("protonBeamDump.material.core"),
            &rotationInShield,
            backSteelPositionInShield,
            beamDumpBack, 0,
            G4Colour::Blue()
            );

    constructExtMonFNAL(beamDumpFront, shieldingRot, parent, CLHEP::HepRotation::IDENTITY, config);

  } // constructProtonBeamDump()

} // namespace mu2e
