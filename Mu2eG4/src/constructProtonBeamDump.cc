// Andrei Gaponenko, 2011

#include "Mu2eG4/inc/constructProtonBeamDump.hh"
#include "Mu2eG4/inc/constructExtMonFNAL.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"

#include <iostream>
#include <cmath>

#include "G4Color.hh"
#include "G4LogicalVolume.hh"
#include "G4Trap.hh"
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4Polycone.hh"
#include "G4ExtrudedSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4TwoVector.hh"
#include "G4SDManager.hh"

#include "G4UniformMagField.hh"
#include "G4Mag_UsualEqRhs.hh"
#include "G4ExactHelixStepper.hh"
//#include "G4NystromRK4.hh"
#include "G4ChordFinder.hh"
#include "G4FieldManager.hh"

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib/exception.h"

#include "GeometryService/inc/GeomHandle.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "Mu2eBuildingGeom/inc/Mu2eBuilding.hh"
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"
#include "GeometryService/inc/WorldG4.hh"
#include "GeometryService/inc/Mu2eEnvelope.hh"

#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/FieldMgr.hh"

#include "Mu2eG4/inc/finishNesting.hh"
#include "G4Helper/inc/VolumeInfo.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  void constructProtonBeamDump(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<ProtonBeamDump> dump;
    GeomHandle<Mu2eBuilding> building;
    GeomHandle<ExtMonFNALBuilding> emfb;
    GeomHandle<Mu2eEnvelope> env;
    GeomHandle<WorldG4> world;

    MaterialFinder materialFinder(config);

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;

    //----------------------------------------------------------------
    // Re-fill a part of the formal "HallAir" with dirt.
    // Use an extruded solid to have a properly angled facet
    // at which the beam dump can be placed.

    std::vector<G4TwoVector> beamDumpDirtOutiline;

    std::copy(building->concreteOuterOutline3().begin(),
              building->concreteOuterOutline3().end(),
              std::back_inserter(beamDumpDirtOutiline));

    std::copy(building->concreteOuterOutline1().begin(),
              building->concreteOuterOutline1().end(),
              std::back_inserter(beamDumpDirtOutiline));

    // points need to be in the clock-wise order, need to reverse:
    std::reverse(beamDumpDirtOutiline.begin(), beamDumpDirtOutiline.end());

    // Add points to complete the outline
    const CLHEP::Hep3Vector mu2eCenterInHall(world->mu2eOriginInWorld() - world->hallFormalCenterInWorld());
    const double hallFormalZminInMu2e  = -world->hallFormalHalfSize()[2] - mu2eCenterInHall.z();
    const double dumpDirtXmin = env->xmin();
    const double dumpDirtXmax = env->xmax();

    if(dumpDirtXmin < beamDumpDirtOutiline.back().x()) {
      beamDumpDirtOutiline.push_back(G4TwoVector(dumpDirtXmin, beamDumpDirtOutiline.back().y()));
    }

    // Two points at the back
    beamDumpDirtOutiline.push_back(G4TwoVector(dumpDirtXmin, hallFormalZminInMu2e));
    beamDumpDirtOutiline.push_back(G4TwoVector(dumpDirtXmax, hallFormalZminInMu2e));

    if(dumpDirtXmax > beamDumpDirtOutiline.front().x()) {
      beamDumpDirtOutiline.push_back(G4TwoVector(dumpDirtXmax, beamDumpDirtOutiline.front().y()));
    }

    // We want to rotate the X'Y' plane of the extruded solid
    // to become XZ plane of Mu2e.   Need to rotate by +90 degrees around X.
    // Because the interfaces below use a backward interpretation of rotations,
    // it's more convenient to work with the inverse matrix
    static const CLHEP::HepRotation beamDumpDirtRotationInv(CLHEP::HepRotationX(-90*CLHEP::degree));

    // To the bottom of the formal hall box, in mu2e coords
    const double dumpDirtYmin =
      world->hallFormalCenterInWorld()[1] - world->hallFormalHalfSize()[1]
      - world->mu2eOriginInWorld()[1]
      ;

    const double dumpDirtYmax = std::max(
                                         emfb->roomInsideYmax() + emfb->roomCeilingThickness() + emfb->dirtOverheadThickness()
                                         ,
                                         dump->frontShieldingCenterInMu2e()[1] + dump->frontShieldingHalfSize()[1]
                                         );


    VolumeInfo beamDumpDirt("ProtonBeamDumpDirt",
                            CLHEP::Hep3Vector(0, (dumpDirtYmax+dumpDirtYmin)/2, 0)
                            - parent.centerInMu2e(),
                            parent.centerInWorld);

    beamDumpDirt.solid = new G4ExtrudedSolid(beamDumpDirt.name, beamDumpDirtOutiline,
                                             (dumpDirtYmax - dumpDirtYmin)/2,
                                             G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

    finishNesting(beamDumpDirt,
                  materialFinder.get("dirt.overburdenMaterialName"),
                  &beamDumpDirtRotationInv,
                  beamDumpDirt.centerInParent,
                  parent.logical,
                  0,
                  config.getBool("protonBeamDump.dirtVisible"),
                  G4Colour(0.9, 0, 0.9), //G4Colour::Magenta(),
                  config.getBool("protonBeamDump.dirtSolid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    //----------------------------------------------------------------
    CLHEP::Hep3Vector frontShieldingPositionInDirt( beamDumpDirtRotationInv * (dump->frontShieldingCenterInMu2e() - beamDumpDirt.centerInMu2e()));

    // finishNesting() uses the backwards interpretation of rotations.
    // We need: (beamDumpDirtRotation.inverse()*dump.coreRotationInMu2e()).inverse()
    static const CLHEP::HepRotation rotationInDirtInv( (beamDumpDirtRotationInv*dump->coreRotationInMu2e()).inverse() );

    // FIXME: from nestBo() here we get VolumeInfo with wrong centerInWorld
    VolumeInfo frontShielding = nestBox("ProtonBeamDumpFrontShielding",
                                        dump->frontShieldingHalfSize(),
                                        materialFinder.get("protonBeamDump.material.shielding"),
                                        &rotationInDirtInv,
                                        frontShieldingPositionInDirt,
                                        beamDumpDirt, 0,
                                        config.getBool("protonBeamDump.shieldingVisible"),
                                        G4Colour::Grey(),
                                        config.getBool("protonBeamDump.shieldingSolid"),
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );

    // FIXME: we should not need to correct the wrong information
    frontShielding.centerInWorld = GeomHandle<WorldG4>()->mu2eOriginInWorld() + dump->frontShieldingCenterInMu2e();


    nestBox("ProtonBeamDumpBackShielding",
            dump->backShieldingHalfSize(),
            materialFinder.get("protonBeamDump.material.shielding"),
            &rotationInDirtInv,
            beamDumpDirtRotationInv * (dump->backShieldingCenterInMu2e() - beamDumpDirt.centerInMu2e()),
            beamDumpDirt, 0,
            config.getBool("protonBeamDump.shieldingVisible"),
            G4Colour::Grey(),
            config.getBool("protonBeamDump.shieldingSolid"),
            forceAuxEdgeVisible,
            placePV,
            doSurfaceCheck
            );

    nestBox("ProtonBeamDumpCore",
            dump->coreHalfSize(),
            materialFinder.get("protonBeamDump.material.core"),
            0,
            dump->coreRotationInMu2e().inverse()*(dump->coreCenterInMu2e() - dump->frontShieldingCenterInMu2e()),
            frontShielding, 0,
            config.getBool("protonBeamDump.coreVisible"),
            G4Colour::Blue(),
            config.getBool("protonBeamDump.coreSolid"),
            forceAuxEdgeVisible,
            placePV,
            doSurfaceCheck
            );

    nestBox("ProtonBeamDumpMouth",
            dump->mouthHalfSize(),
            materialFinder.get("protonBeamDump.material.air"),
            0,
            dump->coreRotationInMu2e().inverse()*(dump->mouthCenterInMu2e() - dump->frontShieldingCenterInMu2e()),
            frontShielding, 0,
            config.getBool("protonBeamDump.mouthVisible"),
            G4Colour::Cyan(),
            config.getBool("protonBeamDump.mouthSolid"),
            forceAuxEdgeVisible,
            placePV,
            doSurfaceCheck
            );

    nestBox("ProtonBeamNeutronCave",
            dump->neutronCaveHalfSize(),
            materialFinder.get("protonBeamDump.material.air"),
            0,
            dump->coreRotationInMu2e().inverse()*(dump->neutronCaveCenterInMu2e() - dump->frontShieldingCenterInMu2e()),
            frontShielding, 0,
            config.getBool("protonBeamDump.neutronCaveVisible"),
            G4Colour::Cyan(),
            config.getBool("protonBeamDump.neutronCaveSolid"),
            forceAuxEdgeVisible,
            placePV,
            doSurfaceCheck
            );

    //----------------------------------------------------------------
    constructExtMonFNAL(frontShielding, dump->coreRotationInMu2e(), beamDumpDirt, beamDumpDirtRotationInv.inverse(), config);

  } // constructProtonBeamDump()

} // namespace mu2e
