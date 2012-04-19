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
#include "GeometryService/inc/WorldG4.hh"

#include "ExtinctionMonitorFNAL/inc/ExtMonFNAL.hh"

#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "Mu2eUtilities/inc/SimpleConfig.hh"
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

  //================================================================
  void constructCollimatorExtMonFNAL(const ProtonBeamDump::CollimatorExtMonFNAL& collimator,
                                     const VolumeInfo& parent,
                                     const CLHEP::Hep3Vector& collimatorCenterInParent,
                                     const CLHEP::HepRotation& collimatorRotationInParent,
                                     const SimpleConfig& config
                                     )
  {
    GeomHandle<ProtonBeamDump> dump;

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
                                      dump->enclosureHalfSize()[0],
                                      dump->enclosureHalfSize()[1],
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
                  config.getBool("extMonFilter."+collimator.name()+".alignmentHole.visible"),
                  G4Colour::Red(),//                  G4Colour::Cyan(),
                  config.getBool("extMonFilter."+collimator.name()+".alignmentHole.solid"),
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
                  config.getBool("extMonFilter."+collimator.name()+".alignmentPlug.visible"),
                  G4Colour(0.4, 0, 0), // G4Colour::Red(),
                  config.getBool("extMonFilter."+collimator.name()+".alignmentPlug.solid"),
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
                  config.getBool("extMonFilter."+collimator.name()+".channel.visible"),
                  G4Colour::Yellow(),
                  config.getBool("extMonFilter."+collimator.name()+".channel.solid"),
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
                  config.getBool("extMonFilter."+collimator.name()+".channel.visible"),
                  G4Colour::Yellow(),
                  config.getBool("extMonFilter."+collimator.name()+".channel.solid"),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

  }

  //================================================================
  void constructFilterMagnet(const ProtonBeamDump& dump,
                             const VolumeInfo& parent,
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
    // We need: (dump.enclosureRotationInMu2e().inverse()*dump.filterMagnetRotationInMu2e()).inverse()
    // == dump.filterMagnetRotationInMu2e().inverse() * dump.enclosureRotationInMu2e()

    CLHEP::HepRotation *magnetRotationInParent =
      reg.add(dump.filterMagnetRotationInMu2e().inverse()*dump.enclosureRotationInMu2e());

    const CLHEP::Hep3Vector nx(1, 0, 0);
    const CLHEP::Hep3Vector ny(0, 1, 0);
    const CLHEP::Hep3Vector nz(0, 0, 1);
    std::cout<<"AG: DEGUB: magnetRotationInParent.inv * Nx = \n"
             <<magnetRotationInParent->inverse()*nx<<"\n"
             <<magnetRotationInParent->inverse()*ny<<"\n"
             <<magnetRotationInParent->inverse()*nz<<"\n"
             <<std::endl;

    const VolumeInfo magnetIron = nestBox("ExtMonFNALFilterMagnetIron",
                                          dump.filterMagnet().outerHalfSize(),
                                          materialFinder.get("extMonFilter.magnet.material"),
                                          magnetRotationInParent,
                                          dump.filterMagnetCenterInEnclosure() - dump.magnetPitCenterInEnclosure(),
                                          parent, 0,
                                          config.getBool("extMonFilter.magnet.iron.visible"),
                                          G4Colour::Magenta(),
                                          config.getBool("extMonFilter.magnet.iron.solid"),
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );

    std::vector<double> apertureHalfSize(3);
    apertureHalfSize[0] = 0.5*dump.filterMagnet().apertureWidth();
    apertureHalfSize[1] = 0.5*dump.filterMagnet().apertureHeight();
    apertureHalfSize[2] = dump.filterMagnet().outerHalfSize()[2];

    VolumeInfo magnetAperture = nestBox("ExtMonFNALFilterMagnetAperture",
                                        apertureHalfSize,
                                        materialFinder.get("hall.insideMaterialName"),
                                        0,
                                        CLHEP::Hep3Vector(0, 0, 0),
                                        magnetIron.logical, 0,
                                        config.getBool("extMonFilter.magnet.aperture.visible"),
                                        G4Colour::Grey(),
                                        config.getBool("extMonFilter.magnet.aperture.solid"),
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );


    //----------------------------------------------------------------
    // Define the field in the magnet

    const CLHEP::Hep3Vector BFieldVector(dump.filterMagnet().fieldStrength()
                                         * (dump.filterMagnetRotationInMu2e()*CLHEP::Hep3Vector(1,0,0)));

    std::cout<<"AG: ExtMonFNAL filter magnet field = "<<BFieldVector<<std::endl;

    G4MagneticField *field = reg.add(new G4UniformMagField(BFieldVector));

    G4Mag_UsualEqRhs *rhs  = reg.add(new G4Mag_UsualEqRhs(field));

    G4MagIntegratorStepper *integrator = reg.add(new G4ExactHelixStepper(rhs));
    //G4MagIntegratorStepper *integrator = reg.add(new G4NystromRK4(rhs));

    const double stepMinimum = config.getDouble("extMonFilter.magnet.stepMinimum", 1.0e-2 * CLHEP::mm /*The default from G4ChordFinder.hh*/);
    G4ChordFinder          *chordFinder = reg.add(new G4ChordFinder(field, stepMinimum, integrator));

    const double deltaOld = chordFinder->GetDeltaChord();
    chordFinder->SetDeltaChord(config.getDouble("extMonFilter.magnet.deltaChord", deltaOld));
    AGDEBUG("chordFinder: using deltaChord = "<<chordFinder->GetDeltaChord()<<" (default = "<<deltaOld<<")");

    G4FieldManager *manager = reg.add(new G4FieldManager(field, chordFinder));

    AGDEBUG("orig: manager epsMin = "<<manager->GetMinimumEpsilonStep()
            <<", epsMax = "<<manager->GetMaximumEpsilonStep()
            <<", deltaOneStep = "<<manager->GetDeltaOneStep()
            );

    manager->SetMinimumEpsilonStep(config.getDouble("extMonFilter.magnet.minEpsilonStep", manager->GetMinimumEpsilonStep()));
    manager->SetMaximumEpsilonStep(config.getDouble("extMonFilter.magnet.maxEpsilonStep", manager->GetMaximumEpsilonStep()));
    manager->SetDeltaOneStep(config.getDouble("extMonFilter.magnet.deltaOneStep", manager->GetDeltaOneStep()));

    AGDEBUG("new:  manager epsMin = "<<manager->GetMinimumEpsilonStep()
            <<", epsMax = "<<manager->GetMaximumEpsilonStep()
            <<", deltaOneStep = "<<manager->GetDeltaOneStep()
            );

    //magnetAperture.logical->SetFieldManager(manager, true);
    magnetIron.logical->SetFieldManager(manager, true);

    //----------------------------------------------------------------
  }

  //================================================================

  void constructProtonBeamDump(const VolumeInfo& parent, const SimpleConfig& config) {

    GeomHandle<ProtonBeamDump> dump;
    GeomHandle<Mu2eBuilding> building;
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

    // Add the last two points to complete the outline
    const CLHEP::Hep3Vector mu2eCenterInHall(world->mu2eOriginInWorld() - world->hallFormalCenterInWorld());
    const double hallFormalZminInMu2e  = -world->hallFormalHalfSize()[2] - mu2eCenterInHall.z();

    beamDumpDirtOutiline.push_back(G4TwoVector(building->hallInsideXmin() - building->hallWallThickness(),
                                               hallFormalZminInMu2e));

    beamDumpDirtOutiline.push_back(G4TwoVector(building->hallInsideXmax() + building->hallWallThickness() ,
                                               hallFormalZminInMu2e));

    // We want to rotate the X'Y' plane of the extruded solid
    // to become XZ plane of Mu2e.   Need to rotate by +90 degrees around X.
    // Because the interfaces below use a backward interpretation of rotations,
    // it's more convenient to work with the inverse matrix
    static CLHEP::HepRotation beamDumpDirtRotationInv(CLHEP::HepRotation::IDENTITY);
    beamDumpDirtRotationInv.rotateX(-90*CLHEP::degree);

    const double dumpDirtYmin = world->dumpDirtFormalYminInMu2e();
    const double dumpDirtYmax = world->dumpDirtFormalYmaxInMu2e();

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
    CLHEP::Hep3Vector enclosurePositionInDirt( beamDumpDirtRotationInv * (dump->enclosureCenterInMu2e() - beamDumpDirt.centerInMu2e()));

    // finishNesting() uses the backwards interpretation of rotations.
    // We need: (beamDumpDirtRotation.inverse()*dump.enclosureRotationInMu2e()).inverse()
    static CLHEP::HepRotation rotationInDirt( (beamDumpDirtRotationInv*dump->enclosureRotationInMu2e()).inverse() );

    const VolumeInfo logicalEnclosure = nestBox("ProtonBeamDumpShielding",
                                                dump->enclosureHalfSize(),
                                                materialFinder.get("protonBeamDump.material.shielding"),
                                                &rotationInDirt,
                                                enclosurePositionInDirt,
                                                beamDumpDirt, 0,
                                                config.getBool("protonBeamDump.logicalEnclosureVisible"),
                                                G4Colour::Grey(),
                                                config.getBool("protonBeamDump.logicalEnclosureSolid"),
                                                forceAuxEdgeVisible,
                                                placePV,
                                                doSurfaceCheck
                                                );

    nestBox("ProtonBeamDumpCore",
            dump->coreHalfSize(),
            materialFinder.get("protonBeamDump.material.core"),
            0,
            dump->coreCenterInEnclosure(),
            logicalEnclosure, 0,
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
            CLHEP::Hep3Vector(0.,
                              dump->coreCenterInEnclosure()[1],
                              dump->enclosureHalfSize()[2] - dump->mouthHalfSize()[2]),
            logicalEnclosure, 0,
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
            CLHEP::Hep3Vector(0.,
                              dump->coreCenterInEnclosure()[1],
                              dump->enclosureHalfSize()[2] - 2*dump->mouthHalfSize()[2] - dump->neutronCaveHalfSize()[2]),
            logicalEnclosure, 0,
            config.getBool("protonBeamDump.neutronCaveVisible"),
            G4Colour::Cyan(),
            config.getBool("protonBeamDump.neutronCaveSolid"),
            forceAuxEdgeVisible,
            placePV,
            doSurfaceCheck
            );

    const VolumeInfo magnetPit  = nestBox("ProtonBeamDumpMagnetPit",
                                          dump->magnetPitHalfSize(),
                                          materialFinder.get("protonBeamDump.material.air"),
                                          0,
                                          dump->magnetPitCenterInEnclosure(),
                                          logicalEnclosure, 0,
                                          config.getBool("protonBeamDump.magnetPitVisible"),
                                          G4Colour::Cyan(),
                                          config.getBool("protonBeamDump.magnetPitSolid"),
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );

    //----------------------------------------------------------------
    // The ExtMon filter channel

    constructCollimatorExtMonFNAL(dump->collimator1(),
                                  logicalEnclosure,
                                  dump->collimator1CenterInEnclosure(),
                                  dump->enclosureRotationInMu2e().inverse()*dump->collimator1RotationInMu2e(),
                                  config);

    constructFilterMagnet(*dump, magnetPit, config);

    constructCollimatorExtMonFNAL(dump->collimator2(),
                                  logicalEnclosure,
                                  dump->collimator2CenterInEnclosure(),
                                  dump->enclosureRotationInMu2e().inverse()*dump->collimator2RotationInMu2e(),
                                  config);

    //----------------------------------------------------------------
    // Add some volumes for visualization purposes.
    //
    // ROOT's OpenGL graphics refuses to display in a nice way the
    // ProtonBeamDumpShielding volume.  On the other hand it turns out
    // the same viewer shows existing VirtualDetectors nicely.
    //
    // Here we outline the dump shielding using very thin boxes (with
    // the same thickness as the VirtualDetector dimension)
    //

    const bool applyVisualizationKludge = config.getBool("protonBeamDump.applyROOTVisualizationKludge", false);
    if(applyVisualizationKludge) {

      const double kludgeHalfThickness = 0.01; // the thickness that works with the current root opengl
      int const nSurfaceCheckPoints = 100000; // for a more thorrow check due to the small thickness

      const bool kludgeIsVisible      = true;
      const bool kludgeIsSolid        = true; //false

      G4Material* vacuumMaterial      = materialFinder.get("toyDS.insideMaterialName");

      // The vertical side walls go inside the dump concrete
      if(1) {

        std::vector<double> hlen(3);
        hlen[0] = kludgeHalfThickness;
        hlen[1] = dump->enclosureHalfSize()[1];
        hlen[2] = dump->enclosureHalfSize()[2];

        VolumeInfo pxInfo = nestBox("BeamDumpShieldingVisKludgePosX",
                                    hlen,
                                    vacuumMaterial,
                                    0,
                                    CLHEP::Hep3Vector(+dump->enclosureHalfSize()[0] - kludgeHalfThickness, 0., 0.),
                                    logicalEnclosure,
                                    0,
                                    kludgeIsVisible,
                                    G4Color::Grey(),
                                    kludgeIsSolid,
                                    forceAuxEdgeVisible,
                                    true,
                                    false /*surface check*/
                                    );

        // the volumes are very thin, a more thorough check is needed
        doSurfaceCheck && pxInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);

        VolumeInfo nxInfo = nestBox("BeamDumpShieldingVisKludgeNegX",
                                    hlen,
                                    vacuumMaterial,
                                    0,
                                    CLHEP::Hep3Vector(-dump->enclosureHalfSize()[0] + kludgeHalfThickness, 0., 0.),
                                    logicalEnclosure,
                                    0,
                                    kludgeIsVisible,
                                    G4Color::Grey(),
                                    kludgeIsSolid,
                                    forceAuxEdgeVisible,
                                    true,
                                    false /*surface check*/
                                    );

        // the volumes are very thin, a more thorough check is needed
        doSurfaceCheck && nxInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
      }

      // The top "visualization plane kludge" can't be put inside because of the magnet pit volume
      // So put both top and bottom planes outside
      //
      // The top surface:
      if(1) {

        std::vector<double> hlen(3);
        hlen[0] = dump->enclosureHalfSize()[0];
        hlen[1] = kludgeHalfThickness;
        hlen[2] = dump->enclosureHalfSize()[2];

        CLHEP::Hep3Vector kludgeCenterInMu2e(dump->enclosureCenterInMu2e()[0],
                                             dump->enclosureCenterInMu2e()[1]+dump->enclosureHalfSize()[1],
                                             dump->enclosureCenterInMu2e()[2]
                                             );

        CLHEP::Hep3Vector kludgeCenterInDirt( beamDumpDirtRotationInv*(kludgeCenterInMu2e - beamDumpDirt.centerInMu2e()) );


        CLHEP::Hep3Vector kludgeOffset(0, 0, kludgeHalfThickness);

        VolumeInfo pyInfo = nestBox("BeamDumpShieldingVisKludgePosY",
                                    hlen,
                                    vacuumMaterial,
                                    &rotationInDirt,
                                    kludgeCenterInDirt + kludgeOffset,
                                    beamDumpDirt, // logicalEnclosure,
                                    0,
                                    kludgeIsVisible,
                                    G4Color::Grey(),
                                    kludgeIsSolid,
                                    forceAuxEdgeVisible,
                                    true,
                                    false /*surface check*/
                                    );

        // the volumes are very thin, a more thorough check is needed
        doSurfaceCheck && pyInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
      }

      // The bottom surface:
      if(1) {

        std::vector<double> hlen(3);
        hlen[0] = dump->enclosureHalfSize()[0];
        hlen[1] = kludgeHalfThickness;
        hlen[2] = dump->enclosureHalfSize()[2];

        CLHEP::Hep3Vector kludgeCenterInMu2e(dump->enclosureCenterInMu2e()[0],
                                             dump->enclosureCenterInMu2e()[1]-dump->enclosureHalfSize()[1],
                                             dump->enclosureCenterInMu2e()[2]
                                             );

        CLHEP::Hep3Vector kludgeCenterInDirt( beamDumpDirtRotationInv*(kludgeCenterInMu2e - beamDumpDirt.centerInMu2e()) );


        CLHEP::Hep3Vector kludgeOffset(0, 0, -kludgeHalfThickness);

        VolumeInfo nyInfo = nestBox("BeamDumpShieldingVisKludgeNegY",
                                    hlen,
                                    vacuumMaterial,
                                    &rotationInDirt,
                                    kludgeCenterInDirt + kludgeOffset,
                                    beamDumpDirt, // logicalEnclosure,
                                    0,
                                    kludgeIsVisible,
                                    G4Color::Grey(),
                                    kludgeIsSolid,
                                    forceAuxEdgeVisible,
                                    true,
                                    false /*surface check*/
                                    );

        // the volumes are very thin, a more thorough check is needed
        doSurfaceCheck && nyInfo.physical->CheckOverlaps(nSurfaceCheckPoints,0.0,true);
      }
    }


    //----------------------------------------------------------------
    if(art::ServiceHandle<GeometryService>()->hasElement<mu2e::ExtMonFNAL::ExtMon>()) {
      constructExtMonFNAL(beamDumpDirt, beamDumpDirtRotationInv.inverse(), config);
    }

  }

}
