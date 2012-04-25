// Andrei Gaponenko, 2012

#include "Mu2eG4/inc/constructExtMonFNAL.hh"

#include <iostream>
#include <cmath>

#include "G4Color.hh"
#include "G4RotationMatrix.hh"
#include "G4LogicalVolume.hh"
#include "G4ExtrudedSolid.hh"
#include "G4Trap.hh"
#include "G4Orb.hh"
#include "G4Box.hh"
#include "G4Polycone.hh"
#include "G4ExtrudedSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4TwoVector.hh"
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
#include "ExtinctionMonitorFNAL/inc/ExtMonFNALBuilding.hh"

#include "G4Helper/inc/VolumeInfo.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/finishNesting.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/nestBox.hh"
#include "Mu2eG4/inc/MaterialFinder.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Mu2eG4/inc/FieldMgr.hh"

// FIXME: should not need that
#include "GeometryService/inc/WorldG4.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

namespace mu2e {

  //================================================================
  void constructCollimatorExtMonFNAL(const ExtMonFNALBuilding::CollimatorExtMonFNAL& collimator,
                                     const VolumeInfo& parent,
                                     const CLHEP::Hep3Vector& collimatorCenterInParent,
                                     const CLHEP::HepRotation& collimatorRotationInParent,
                                     double cutboxdx, double cutboxdy,
                                     const SimpleConfig& config
                                     )
  {
    GeomHandle<ExtMonFNALBuilding> emfb;

    MaterialFinder materialFinder(config);
    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;

    // The G4 interface we use through finishNesting() applies
    // the backwards interpretation of rotations.  Be consistent
    // and use the "backwards" interface of boolean solids.
    CLHEP::HepRotation *colrot = reg.add(collimatorRotationInParent.inverse());

    // Make sure the solids definining the collimator hole etc are sufficiently long to completely
    // exit the concrete on both ends.

    using std::abs;
    const double boxdz = 0.5*collimator.horizontalLength();
    const double dr =  *std::max_element(collimator.alignmentPlugRadius().begin(), collimator.alignmentPlugRadius().end());
    const double cylHalfLength = std::sqrt(std::pow(boxdz, 2) +
                                           std::pow(boxdz * abs(tan(collimator.angleV())/cos(collimator.angleH())) + dr/abs(cos(collimator.angleV())), 2) +
                                           std::pow(boxdz * abs(tan(collimator.angleH())) + dr/abs(cos(collimator.angleH())), 2)
                                           );

    double zPlane[] = {-cylHalfLength, -0.5*collimator.radiusTransitiondZ(), +0.5*collimator.radiusTransitiondZ(), +cylHalfLength };
    double rzero[] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};


    G4Box* cutbox = reg.add(new G4Box(collimator.name()+"cutbox",
                                      cutboxdx,
                                      cutboxdy,
                                      0.5*collimator.horizontalLength())
                            );

    //----------------------------------------------------------------
    // Alignment hole

    double rAlignmentHole[] = {
      collimator.alignmentPlugRadius()[1] + collimator.alignmentHoleRClearance()[1],
      collimator.alignmentPlugRadius()[1] + collimator.alignmentHoleRClearance()[1],
      collimator.alignmentPlugRadius()[0] + collimator.alignmentHoleRClearance()[0],
      collimator.alignmentPlugRadius()[0] + collimator.alignmentHoleRClearance()[0]
    };

    G4Polycone *holeCylinder = reg.add(new G4Polycone(collimator.name()+"holecomponent", 0, 2*M_PI,
                                                      sizeof(zPlane)/sizeof(zPlane[0]),
                                                      zPlane,
                                                      rzero,
                                                      rAlignmentHole
                                                      )
                                       );


    VolumeInfo alignmentHole(collimator.name()+"AlignmentHole",
                             collimatorCenterInParent,
                             parent.centerInWorld);

    alignmentHole.solid = new G4IntersectionSolid(alignmentHole.name,
                                                  cutbox,
                                                  holeCylinder,
                                                  colrot,
                                                  CLHEP::Hep3Vector(0,0,0)
                                                  );

    finishNesting(alignmentHole,
                  materialFinder.get("hall.insideMaterialName"),
                  0,
                  alignmentHole.centerInParent,
                  parent.logical,
                  0,
                  config.getBool("extMonFNAL."+collimator.name()+".alignmentHole.visible"),
                  G4Colour::Red(),//                  G4Colour::Cyan(),
                  config.getBool("extMonFNAL."+collimator.name()+".alignmentHole.solid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    //----------------------------------------------------------------
    // Alignment plug

    double rAlignmentPlug[] = {
      collimator.alignmentPlugRadius()[1],
      collimator.alignmentPlugRadius()[1],
      collimator.alignmentPlugRadius()[0],
      collimator.alignmentPlugRadius()[0]
    };

    G4Polycone *plugCylinder = reg.add(new G4Polycone(collimator.name()+"plugcomponent", 0, 2*M_PI,
                                                      sizeof(zPlane)/sizeof(zPlane[0]),
                                                      zPlane,
                                                      rzero,
                                                      rAlignmentPlug
                                                      )
                                       );

    VolumeInfo alignmentPlug(collimator.name()+"AlignmentPlug",
                             CLHEP::Hep3Vector(0,0,0),
                             alignmentHole.centerInWorld);

    alignmentPlug.solid = new G4IntersectionSolid(alignmentPlug.name,
                                                  cutbox,
                                                  plugCylinder,
                                                  colrot,
                                                  CLHEP::Hep3Vector(0,0,0)
                                                  );

    finishNesting(alignmentPlug,
                  materialFinder.get("protonBeamDump.material.shielding"),
                  0,
                  alignmentPlug.centerInParent,
                  alignmentHole.logical,
                  0,
                  config.getBool("extMonFNAL."+collimator.name()+".alignmentPlug.visible"),
                  G4Colour(0.4, 0, 0), // G4Colour::Red(),
                  config.getBool("extMonFNAL."+collimator.name()+".alignmentPlug.solid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    //----------------------------------------------------------------
    // Collimator channel: the upstream half

    G4Box *upbox = reg.add(new G4Box(collimator.name()+"upbox",
                                     0.5*collimator.channelWidth()[0],
                                     0.5*collimator.channelHeight()[0],
                                     0.5*cylHalfLength
                                     )
                           );

    VolumeInfo channelUp(collimator.name()+"ChannelUp",
                         CLHEP::Hep3Vector(0, 0, 0),
                         alignmentPlug.centerInWorld);

    channelUp.solid = new G4IntersectionSolid(channelUp.name,
                                              cutbox,
                                              upbox,
                                              colrot,
                                              colrot->inverse()*CLHEP::Hep3Vector(0,0,0.5*cylHalfLength)
                                              );

    finishNesting(channelUp,
                  materialFinder.get("hall.insideMaterialName"),
                  0,
                  channelUp.centerInParent,
                  alignmentPlug.logical,
                  0,
                  config.getBool("extMonFNAL."+collimator.name()+".channel.visible"),
                  G4Colour::Yellow(),
                  config.getBool("extMonFNAL."+collimator.name()+".channel.solid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    //----------------------------------------------------------------
    // Collimator channel: the downstream half

    G4Box *downbox = reg.add(new G4Box(collimator.name()+"downbox",
                                       0.5*collimator.channelWidth()[1],
                                       0.5*collimator.channelHeight()[1],
                                       0.5*cylHalfLength
                                       )
                             );

    VolumeInfo channelDown(collimator.name()+"ChannelDown",
                           CLHEP::Hep3Vector(0, 0, 0),
                           alignmentPlug.centerInWorld);

    channelDown.solid = new G4IntersectionSolid(channelDown.name,
                                                cutbox,
                                                downbox,
                                                colrot,
                                                colrot->inverse()*CLHEP::Hep3Vector(0, 0, -0.5*cylHalfLength)
                                                );

    finishNesting(channelDown,
                  materialFinder.get("hall.insideMaterialName"),
                  0,
                  channelDown.centerInParent,
                  alignmentPlug.logical,
                  0,
                  config.getBool("extMonFNAL."+collimator.name()+".channel.visible"),
                  G4Colour::Yellow(),
                  config.getBool("extMonFNAL."+collimator.name()+".channel.solid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

  }

  //================================================================
  void constructFilterMagnet(const ExtMonFNALBuilding& emfb,
                             const VolumeInfo& parent,
                             const CLHEP::HepRotation& parentRotationInMu2e,
                             const SimpleConfig& config
                             )
  {
    MaterialFinder materialFinder(config);

    const bool forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible");
    const bool doSurfaceCheck      = config.getBool("g4.doSurfaceCheck");
    const bool placePV             = true;

    AntiLeakRegistry& reg = art::ServiceHandle<G4Helper>()->antiLeakRegistry();

    //----------------------------------------------------------------
    // finishNesting() uses the backwards interpretation of rotations
    CLHEP::HepRotation *magnetRotationInParentInv =
      // (parentRotationInMu2e.inverse()*magnetRotationInMu2e).inverse()
      reg.add(emfb.filterMagnetRotationInMu2e().inverse()*parentRotationInMu2e);

    const CLHEP::Hep3Vector nx(1, 0, 0);
    const CLHEP::Hep3Vector ny(0, 1, 0);
    const CLHEP::Hep3Vector nz(0, 0, 1);
    std::cout<<"AG: DEGUB: magnetRotationInParent.inv * Nx = \n"
             <<*magnetRotationInParentInv*nx<<"\n"
             <<*magnetRotationInParentInv*ny<<"\n"
             <<*magnetRotationInParentInv*nz<<"\n"
             <<std::endl;


    AGDEBUG("emfb.filterMagnetCenterInMu2e() = "<<emfb.filterMagnetCenterInMu2e());
    AGDEBUG("magnet parent.centerInMu2e() = "<<parent.centerInMu2e()<<", in world = "<<parent.centerInWorld<<", in parent = "<<parent.centerInParent);
    AGDEBUG("magnet center in parent = "<<parentRotationInMu2e.inverse()*(emfb.filterMagnetCenterInMu2e() - parent.centerInMu2e()));

    const VolumeInfo magnetIron = nestBox("ExtMonFNALFilterMagnetIron",
                                          emfb.filterMagnet().outerHalfSize(),
                                          materialFinder.get("extMonFNAL.magnet.material"),
                                          magnetRotationInParentInv,
                                          parentRotationInMu2e.inverse()*(emfb.filterMagnetCenterInMu2e() - parent.centerInMu2e()),
                                          parent, 0,
                                          config.getBool("extMonFNAL.magnet.iron.visible"),
                                          G4Colour::Magenta(),
                                          config.getBool("extMonFNAL.magnet.iron.solid"),
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );

    std::vector<double> apertureHalfSize(3);
    apertureHalfSize[0] = 0.5*emfb.filterMagnet().apertureWidth();
    apertureHalfSize[1] = 0.5*emfb.filterMagnet().apertureHeight();
    apertureHalfSize[2] = emfb.filterMagnet().outerHalfSize()[2];

    VolumeInfo magnetAperture = nestBox("ExtMonFNALFilterMagnetAperture",
                                        apertureHalfSize,
                                        materialFinder.get("hall.insideMaterialName"),
                                        0,
                                        CLHEP::Hep3Vector(0, 0, 0),
                                        magnetIron.logical, 0,
                                        config.getBool("extMonFNAL.magnet.aperture.visible"),
                                        G4Colour::Grey(),
                                        config.getBool("extMonFNAL.magnet.aperture.solid"),
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );


    //----------------------------------------------------------------
    // Define the field in the magnet

    const CLHEP::Hep3Vector BFieldVector(emfb.filterMagnet().fieldStrength()
                                         * (emfb.filterMagnetRotationInMu2e()*CLHEP::Hep3Vector(1,0,0)));

    std::cout<<"AG: ExtMonFNAL filter magnet field = "<<BFieldVector<<std::endl;

    G4MagneticField *field = reg.add(new G4UniformMagField(BFieldVector));

    G4Mag_UsualEqRhs *rhs  = reg.add(new G4Mag_UsualEqRhs(field));

    G4MagIntegratorStepper *integrator = reg.add(new G4ExactHelixStepper(rhs));
    //G4MagIntegratorStepper *integrator = reg.add(new G4NystromRK4(rhs));

    const double stepMinimum = config.getDouble("extMonFNAL.magnet.stepMinimum", 1.0e-2 * CLHEP::mm /*The default from G4ChordFinder.hh*/);
    G4ChordFinder          *chordFinder = reg.add(new G4ChordFinder(field, stepMinimum, integrator));

    const double deltaOld = chordFinder->GetDeltaChord();
    chordFinder->SetDeltaChord(config.getDouble("extMonFNAL.magnet.deltaChord", deltaOld));
    AGDEBUG("chordFinder: using deltaChord = "<<chordFinder->GetDeltaChord()<<" (default = "<<deltaOld<<")");

    G4FieldManager *manager = reg.add(new G4FieldManager(field, chordFinder));

    AGDEBUG("orig: manager epsMin = "<<manager->GetMinimumEpsilonStep()
            <<", epsMax = "<<manager->GetMaximumEpsilonStep()
            <<", deltaOneStep = "<<manager->GetDeltaOneStep()
            );

    manager->SetMinimumEpsilonStep(config.getDouble("extMonFNAL.magnet.minEpsilonStep", manager->GetMinimumEpsilonStep()));
    manager->SetMaximumEpsilonStep(config.getDouble("extMonFNAL.magnet.maxEpsilonStep", manager->GetMaximumEpsilonStep()));
    manager->SetDeltaOneStep(config.getDouble("extMonFNAL.magnet.deltaOneStep", manager->GetDeltaOneStep()));

    AGDEBUG("new:  manager epsMin = "<<manager->GetMinimumEpsilonStep()
            <<", epsMax = "<<manager->GetMaximumEpsilonStep()
            <<", deltaOneStep = "<<manager->GetDeltaOneStep()
            );

    //magnetAperture.logical->SetFieldManager(manager, true);
    magnetIron.logical->SetFieldManager(manager, true);
  }

  //================================================================
  VolumeInfo constructExtMonFNALBuilding(const VolumeInfo& collimator1Parent,
                                         const CLHEP::HepRotation& collimator1ParentRotationInMu2e,
                                         const VolumeInfo& mainParent,
                                         const CLHEP::HepRotation& mainParentRotationInMu2e,
                                         const SimpleConfig& config)
  {
    bool const forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    bool const doSurfaceCheck      = config.getBool("g4.doSurfaceCheck",false);
    bool const placePV             = true;

    MaterialFinder materialFinder(config);

    GeomHandle<ProtonBeamDump> dump;
    GeomHandle<ExtMonFNALBuilding> emfb;

    static const CLHEP::HepRotation roomRotationInParentInv((mainParentRotationInMu2e.inverse() * emfb->roomRotationInMu2e()).inverse());

    // Test
    if(false) {
      G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));
      const VolumeInfo& hall = _helper->locateVolInfo("HallAir");
      VolumeInfo test("emfroomtest", emfb->roomRefPointInMu2e() - hall.centerInMu2e(), hall.centerInWorld);
      test.solid = new G4Orb(test.name, 3000.);

      finishNesting(test,
                    materialFinder.get("extMonFNAL.room.materialName"),
                    0,
                    test.centerInParent,
                    hall.logical,
                    0,
                    true/*visible*/,
                    G4Colour::Magenta(),
                    true/*solid*/,
                    forceAuxEdgeVisible,
                    placePV,
                    false
                    );
    }

    const CLHEP::Hep3Vector roomRefPointInParent(mainParentRotationInMu2e.inverse()*(emfb->roomRefPointInMu2e() - mainParent.centerInMu2e()));

    VolumeInfo roomAir("ExtMonFNALRoomAir", roomRefPointInParent, mainParent.centerInWorld);
    roomAir.solid = new G4ExtrudedSolid(roomAir.name,
                                        emfb->roomInsideOutline(),
                                        0.5*emfb->roomInsideFullHeight(),
                                        G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

    // FIXME: we should not need to correct the wrong information
    roomAir.centerInWorld = GeomHandle<WorldG4>()->mu2eOriginInWorld() + emfb->roomRefPointInMu2e();

    finishNesting(roomAir,
                  materialFinder.get("extMonFNAL.room.materialName"),
                  &roomRotationInParentInv,
                  roomAir.centerInParent,
                  mainParent.logical,
                  0,
                  config.getBool("extMonFNAL.room.visible"),
                  G4Colour::Cyan(),
                  config.getBool("extMonFNAL.room.solid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    //----------------------------------------------------------------
    // Collimator2 shielding

    static const CLHEP::HepRotation coll2ShieldingRotationInRoomInv
      (emfb->coll2ShieldingRotationInMu2e().inverse() * emfb->roomRotationInMu2e());

    const CLHEP::Hep3Vector coll2ShieldingCenterInRoom
      (emfb->roomRotationInMu2e().inverse()*(emfb->coll2ShieldingCenterInMu2e() - emfb->roomRefPointInMu2e()));

    VolumeInfo coll2Shielding = nestBox("ExtMonFNALColl2Shielding",
                                        emfb->coll2ShieldingHalfSize(),
                                        materialFinder.get("extMonFNAL.room.wallMaterialName"),
                                        &coll2ShieldingRotationInRoomInv,
                                        coll2ShieldingCenterInRoom,
                                        roomAir,
                                        0,
                                        config.getBool("extMonFNAL.collimator2.shielding.visible"),
                                        G4Colour::Grey(),
                                        config.getBool("extMonFNAL.collimator2.shielding.solid"),
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck);

    // FIXME: we should not need to correct the wrong information
    coll2Shielding.centerInWorld = GeomHandle<WorldG4>()->mu2eOriginInWorld() + emfb->coll2ShieldingCenterInMu2e();

    //----------------------------------------------------------------
    // The filter channel

    AGDEBUG("emfb->collimator1CenterInMu2e() = "<<emfb->collimator1CenterInMu2e());
    AGDEBUG("collimator1Parent.centerInMu2e() = "<<collimator1Parent.centerInMu2e());

    constructCollimatorExtMonFNAL(emfb->collimator1(),
                                  collimator1Parent,

                                  collimator1ParentRotationInMu2e.inverse()*(emfb->collimator1CenterInMu2e()-collimator1Parent.centerInMu2e()),

                                  collimator1ParentRotationInMu2e.inverse()*emfb->collimator1RotationInMu2e(),

                                  dump->frontShieldingHalfSize()[0],
                                  dump->frontShieldingHalfSize()[1],
                                  config);

    constructFilterMagnet(*emfb, roomAir, emfb->roomRotationInMu2e(), config);

    constructCollimatorExtMonFNAL(emfb->collimator2(),
                                  coll2Shielding,
                                  emfb->coll2ShieldingRotationInMu2e().inverse()*(emfb->collimator2CenterInMu2e() - emfb->coll2ShieldingCenterInMu2e()),
                                  emfb->coll2ShieldingRotationInMu2e().inverse()*emfb->collimator2RotationInMu2e(),
                                  emfb->coll2ShieldingHalfSize()[0],
                                  emfb->coll2ShieldingHalfSize()[1],
                                  config);

    // Test
    if(false) {
      G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));
      const VolumeInfo& hall = _helper->locateVolInfo("HallAir");
      VolumeInfo test("emfMagnettest", emfb->filterMagnetCenterInMu2e() - hall.centerInMu2e(), hall.centerInWorld);
      test.solid = new G4Orb(test.name, 500.);

      finishNesting(test,
                    materialFinder.get("extMonFNAL.room.materialName"),
                    0,
                    test.centerInParent,
                    hall.logical,
                    0,
                    true/*visible*/,
                    G4Colour::Blue(),
                    true /*config.getBool("extMonFNAL.roomSolid")*/,
                    forceAuxEdgeVisible,
                    placePV,
                    false
                    );

    }

    //----------------------------------------------------------------
    return roomAir;

  } // constructExtMonFNALBuilding()
} // namespace mu2e
