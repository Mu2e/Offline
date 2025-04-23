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
    geomOptions->loadEntry( config, "ProtonBeamDumpConcrete" , "protonBeamDump.concrete");
    geomOptions->loadEntry( config, "ProtonBeamDumpDirt"     , "protonBeamDump.dirt"       );
    geomOptions->loadEntry( config, "ProtonBeamDumpCore"     , "protonBeamDump.core"       );
    geomOptions->loadEntry( config, "ProtonBeamDumpMouth"    , "protonBeamDump.mouth"      );
    geomOptions->loadEntry( config, "ProtonBeamNeutronCave"  , "protonBeamDump.neutronCave");

    // Use the same rotation as for other concrete volumes in constructHall()
    CLHEP::HepRotation *pshieldingRot = reg.add(new CLHEP::HepRotation());
    CLHEP::HepRotation& shieldingRot = *pshieldingRot;
    pshieldingRot->rotateX( 90*CLHEP::degree);
    pshieldingRot->rotateZ( 90*CLHEP::degree);

    //--------------------------------------------------------------------
    // Subtraction Cylinder

    CLHEP::Hep3Vector subCylOffsetInParent = shieldingRot *(emfb->filter().collimator1().centerInMu2e()
                                                            - dump->dumpConcreteCenterInMu2e());

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
    //Beam Dump Concrete

    static CLHEP::HepRotation rotationInShield( (shieldingRot*dump->coreRotationInMu2e()).inverse() );

    CLHEP::Hep3Vector extMonSubtractionPositionInShield( shieldingRot * (dump->extMonSubtractionCenterInMu2e()
                                                                         - dump->dumpConcreteCenterInMu2e()));

    CLHEP::Hep3Vector mouthPositionInShield( shieldingRot * (dump->mouthCenterInMu2e()
                                                             - dump->dumpConcreteCenterInMu2e()));

    CLHEP::Hep3Vector cavePositionInShield( shieldingRot * (dump->neutronCaveCenterInMu2e()
                                                            - dump->dumpConcreteCenterInMu2e()));

    VolumeInfo beamDumpConcrete("ProtonBeamDumpConcrete",
                                dump->dumpConcreteCenterInMu2e() - parent.centerInMu2e(),
                                parent.centerInWorld
                                );

    G4ExtrudedSolid* beamDumpConcreteExtrusion  = new G4ExtrudedSolid("ProtonBeamDumpConcreteExtrusion",
                                                                      dump->dumpConcreteOutline(),
                                                                      dump->dumpConcreteHalfHeight(),
                                                                      G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.
                                                                      );

    G4ExtrudedSolid* extMonCutout  = new G4ExtrudedSolid("ProtonBeamDumpExtMonCutout",
                                                         dump->extMonSubtractionOutline(),
                                                         dump->extMonSubtractionHalfHeight(),
                                                         G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.
                                                         );

    G4SubtractionSolid* beamDumpConcreteSansExtMon = new G4SubtractionSolid("beamDumpConcreteSansExtMon",
                                                                            beamDumpConcreteExtrusion,
                                                                            extMonCutout,
                                                                            0,
                                                                            extMonSubtractionPositionInShield
                                                                            );

    G4SubtractionSolid* beamDumpConcreteCylSubtraction = new G4SubtractionSolid("beamDumpConcreteCylSubtraction",
                                                                                beamDumpConcreteSansExtMon,
                                                                                subCylinder,
                                                                                subCylinderRotation,
                                                                                subCylOffsetInParent
                                                                                );

    G4Box* ProtonBeamNeutronCave = new G4Box("ProtonBeamNeutronCave",
                                             dump->neutronCaveHalfSize()[0],
                                             dump->neutronCaveHalfSize()[1],
                                             dump->neutronCaveHalfSize()[2]
                                             );

    G4SubtractionSolid* beamDumpConcreteCaveSubtraction = new G4SubtractionSolid("beamDumpConcreteCaveSubtraction",
                                                                                 beamDumpConcreteCylSubtraction,
                                                                                 ProtonBeamNeutronCave,
                                                                                 &rotationInShield,
                                                                                 cavePositionInShield
                                                                                 );

    G4Box* ProtonBeamDumpMouth = new G4Box("ProtonBeamDumpMouth",
                                           dump->mouthHalfSize()[0],
                                           dump->mouthHalfSize()[1],
                                           1.1*dump->mouthHalfSize()[2] //scale up to avoid coinciding boundaries
                                           );

    beamDumpConcrete.solid = new G4SubtractionSolid(beamDumpConcrete.name,
                                                    beamDumpConcreteCaveSubtraction,
                                                    ProtonBeamDumpMouth,
                                                    &rotationInShield,
                                                    mouthPositionInShield
                                                    );

    finishNesting(beamDumpConcrete,
                  materialFinder.get("protonBeamDump.material.shielding"),
                  pshieldingRot,
                  beamDumpConcrete.centerInParent,
                  parent.logical,
                  0,
                  G4Colour::Red() ,
                  "ProtonBeamDumpConcrete"
                  );

    //--------------------------------------------------------------------------------
    CLHEP::Hep3Vector coreAirPositionInShield( shieldingRot * (dump->coreAirCenterInMu2e() - beamDumpConcrete.centerInMu2e()));

    VolumeInfo  beamDumpCoreAir = nestBox("ProtonBeamDumpCoreAir",
                                          dump->coreAirHalfSize(),
                                          materialFinder.get("protonBeamDump.material.air"),
                                          &rotationInShield,
                                          coreAirPositionInShield,
                                          beamDumpConcrete, 0,
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

    CLHEP::Hep3Vector topSteelFlatPositionInShield( shieldingRot * (dump->topSteelFlatCenterInMu2e() - beamDumpConcrete.centerInMu2e()));
    nestBox("ProtonBeamDumpTopSteelFlat",
            dump->topSteelFlatHalfSize(),
            materialFinder.get("protonBeamDump.material.core"),
            &rotationInShield,
            topSteelFlatPositionInShield,
            beamDumpConcrete, 0,
            G4Colour::Blue()
            );

    //----------------------------------------------------------------
    // The scalloped steel

    G4Tubs* scallopCutout = new G4Tubs("BeamDumpScallopCutout",
                                       0.*CLHEP::mm,

                                       emfb->filter().collimator1().shotLinerOuterRadius()
                                       + dump->scallopDistanceToCollimator(),

                                       0.5*emfb->filter().collimator1().length(),
                                       0,
                                       CLHEP::twopi);

    G4Box* scallopedSteelInitialBox  = new G4Box("BeamDumpScallopedSteelInitialBox",
                                                 dump->topSteelScallopedHalfSize()[0],
                                                 dump->topSteelScallopedHalfSize()[1],
                                                 dump->topSteelScallopedHalfSize()[2]
                                                 );

    // The calculation below looks confusing because of the naming
    // correspondance to active vs passive rotations.
    // To get the cutout rotation C in steel S, we need Cinv*S.  But
    // then we want the inverse of it for the G4 interface, leading to
    // (Cinv*S)inv = Sinv * C.  The first two operators below are parts of
    // Sinv.  Then the collimator rotation value has the opposite meaning
    // compared to S, so it also gets the inverse() call below.
    CLHEP::HepRotation *scallopCutoutRotation = reg.add(new CLHEP::HepRotation());
    *scallopCutoutRotation *= shieldingRot.inverse();
    *scallopCutoutRotation *= rotationInShield.inverse();
    *scallopCutoutRotation *= emfb->filter().collimator1().rotationInMu2e().inverse();

    CLHEP::Hep3Vector topSteelScallopedPositionInShield( shieldingRot
                                                         *(dump->topSteelScallopedCenterInMu2e()
                                                           - beamDumpConcrete.centerInMu2e())
                                                         );

    VolumeInfo scallopedSteel("ProtonBeamDumpTopSteelScalloped",
                              topSteelScallopedPositionInShield,
                              parent.centerInWorld);

    CLHEP::Hep3Vector scallopCutoutPositionInSteel( rotationInShield*shieldingRot
                                                    *(emfb->filter().collimator1().centerInMu2e()
                                                      - dump->topSteelScallopedCenterInMu2e()
                                                      )
                                                    );

    scallopedSteel.solid = new G4SubtractionSolid(scallopedSteel.name,
                                                  scallopedSteelInitialBox,
                                                  scallopCutout,
                                                  scallopCutoutRotation,
                                                  scallopCutoutPositionInSteel
                                                  );

    finishNesting(scallopedSteel,
                  materialFinder.get("protonBeamDump.material.core"),
                  &rotationInShield,
                  topSteelScallopedPositionInShield,
                  beamDumpConcrete.logical,
                  0,
                  G4Colour::Blue(),
                  "ProtonBeamDumpTopSteel"
                  );

    //----------------------------------------------------------------
    constructExtMonFNAL(beamDumpConcrete, shieldingRot, parent, CLHEP::HepRotation::IDENTITY, config);

  } // constructProtonBeamDump()

} // namespace mu2e
