// Andrei Gaponenko, 2011

#include "Offline/Mu2eG4/inc/constructProtonBeamDump.hh"
#include "Offline/Mu2eG4/inc/constructExtMonFNAL.hh"
#include "Offline/Mu2eG4Helper/inc/AntiLeakRegistry.hh"

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
#include "Geant4/G4SubtractionSolid.hh"
#include "Geant4/G4TwoVector.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4Tubs.hh"

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

#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"
#include "Offline/GeometryService/inc/WorldG4.hh"
#include "Offline/GeometryService/inc/Mu2eEnvelope.hh"

#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/FieldMgr.hh"

#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  void constructProtonBeamDump(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<ProtonBeamDump> dump;
    GeomHandle<ExtMonFNALBuilding> emfb;
    GeomHandle<Mu2eEnvelope> env;
    GeomHandle<WorldG4> world;

    MaterialFinder materialFinder(config);
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "ProtonBeamDumpFront"    , "protonBeamDump.front");
    geomOptions->loadEntry( config, "ProtonBeamDumpBack"     , "protonBeamDump.back");
    geomOptions->loadEntry( config, "ProtonBeamDumpDirt"     , "protonBeamDump.dirt"       );
    geomOptions->loadEntry( config, "ProtonBeamDumpShielding", "protonBeamDump.shielding"  );
    geomOptions->loadEntry( config, "ProtonBeamDumpCore"     , "protonBeamDump.core"       );
    geomOptions->loadEntry( config, "ProtonBeamDumpMouth"    , "protonBeamDump.mouth"      );
    geomOptions->loadEntry( config, "ProtonBeamNeutronCave"  , "protonBeamDump.neutronCave");

    CLHEP::HepRotation *pshieldingRot = reg.add(new CLHEP::HepRotation());
    CLHEP::HepRotation& shieldingRot = *pshieldingRot;
    pshieldingRot->rotateX( 90*CLHEP::degree);
    pshieldingRot->rotateZ( 90*CLHEP::degree);

    //--------------------------------------------------------------------
    // Subtraction Cylinder

    CLHEP::Hep3Vector subCylOffsetInParent = shieldingRot *(emfb->filter().collimator1().centerInMu2e()
                                                            - CLHEP::Hep3Vector(0, dump->frontShieldingCenterInMu2e()[1], 0));

    //Create a Rotation Matrix
    CLHEP::HepRotation *subCylinderRotation = reg.add(new CLHEP::HepRotation());
    subCylinderRotation->rotateY(90*CLHEP::degree);
    subCylinderRotation->rotateX(emfb->filter().collimator1().angleH_inBeamDump()-dump->coreRotY());
    subCylinderRotation->rotateY(-emfb->filter().collimator1().angleV());

    G4Tubs* subCylinder = new G4Tubs("ExtMonFNALCollimator1Hole",
                                      0.*CLHEP::mm,
                                      emfb->filter().collimator1().shotLinerOuterRadius(),
                                     // factor of 2 to help the hole completely pierce the dump
                                     // even when the collimator center is off
                                     2. * 0.5*emfb->filter().collimator1().length(),
                                      0,
                                      CLHEP::twopi);

    //--------------------------------------------------------------------------------
    //Beam Dump Front

    static CLHEP::HepRotation rotationInShield( (shieldingRot*dump->coreRotationInMu2e()).inverse() );

    CLHEP::Hep3Vector mouthPositionInShield( shieldingRot * (dump->mouthCenterInMu2e()
                                             - (CLHEP::Hep3Vector(0, dump->frontShieldingCenterInMu2e()[1], 0))));

    CLHEP::Hep3Vector cavePositionInShield( shieldingRot * (dump->neutronCaveCenterInMu2e()
                                            - (CLHEP::Hep3Vector(0, dump->frontShieldingCenterInMu2e()[1], 0))));

    VolumeInfo beamDumpFront("ProtonBeamDumpFront",
                             CLHEP::Hep3Vector(0, dump->frontShieldingCenterInMu2e()[1], 0)
                             - parent.centerInMu2e(),
                             parent.centerInWorld
                             );

   G4ExtrudedSolid* beamDumpFrontExtrusion  = new G4ExtrudedSolid("ProtonBeamDumpFrontExtrusion",
                                                                  dump->frontShieldingOutline(),
                                                                  dump->frontShieldingHalfSize()[1],
                                                                  G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.
                                                                  );

   G4SubtractionSolid* beamDumpFrontCylSubtraction = new G4SubtractionSolid("beamDumpFrontCylSubtraction",
                                                                            beamDumpFrontExtrusion,
                                                                            subCylinder,
                                                                            subCylinderRotation,
                                                                            subCylOffsetInParent
                                                                            );

   G4Box* ProtonBeamNeutronCave = new G4Box("ProtonBeamNeutronCave",
                                            dump->neutronCaveHalfSize()[0],
                                            dump->neutronCaveHalfSize()[1],
                                            dump->neutronCaveHalfSize()[2]
                                            );

   G4SubtractionSolid* beamDumpFrontCaveSubtraction = new G4SubtractionSolid("beamDumpFrontCaveSubtraction",
                                                                             beamDumpFrontCylSubtraction,
                                                                             ProtonBeamNeutronCave,
                                                                             &rotationInShield,
                                                                             cavePositionInShield
                                                                             );

   G4Box* ProtonBeamDumpMouth = new G4Box("ProtonBeamDumpMouth",
                                          dump->mouthHalfSize()[0],
                                          dump->mouthHalfSize()[1],
                                          1.1*dump->mouthHalfSize()[2] //scale up to avoid coinciding boundaries
                                          );

   beamDumpFront.solid = new G4SubtractionSolid(beamDumpFront.name,
                                                beamDumpFrontCaveSubtraction,
                                                ProtonBeamDumpMouth,
                                                &rotationInShield,
                                                mouthPositionInShield
                                                );

    finishNesting(beamDumpFront,
                  materialFinder.get("protonBeamDump.material.shielding"),
                  pshieldingRot,
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

    //--------------------------------------------------------------------------------
    //Beam Dump Back

    CLHEP::Hep3Vector coreAirPositionInShield( shieldingRot * (dump->coreAirCenterInMu2e() - beamDumpFront.centerInMu2e()));
    VolumeInfo beamDumpBack("ProtonBeamDumpBack",
                            CLHEP::Hep3Vector(0, dump->backShieldingCenterInMu2e()[1], 0)
                            - parent.centerInMu2e(),
                            parent.centerInWorld);

    beamDumpBack.solid = new G4ExtrudedSolid(beamDumpBack.name, dump->backShieldingOutline(),
                                             dump->backShieldingHalfSize()[1],
                                             G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

    finishNesting(beamDumpBack,
                  materialFinder.get("protonBeamDump.material.shielding"),
                  pshieldingRot,
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
