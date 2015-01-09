// Andrei Gaponenko, 2012

#include "Mu2eG4/inc/constructExtMonFNAL.hh"

#include <iostream>
#include <iterator>
#include <algorithm>
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
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"

#include "G4Helper/inc/VolumeInfo.hh"
#include "G4Helper/inc/G4Helper.hh"
#include "G4Helper/inc/AntiLeakRegistry.hh"
#include "ConfigTools/inc/SimpleConfig.hh"
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

    //----------------------------------------------------------------
    // Alignment hole

    double rAlignmentHole[] = {
      collimator.alignmentHoleRadius()[1],
      collimator.alignmentHoleRadius()[1],
      collimator.alignmentHoleRadius()[0],
      collimator.alignmentHoleRadius()[0]
    };

    G4Polycone *holeCylinder = new G4Polycone(collimator.name()+"holecomponent", 0, 2*M_PI,
                                              sizeof(zPlane)/sizeof(zPlane[0]),
                                              zPlane,
                                              rzero,
                                              rAlignmentHole);


    VolumeInfo alignmentHole(collimator.name()+"AlignmentHole",
                             collimatorCenterInParent,
                             parent.centerInWorld);

    alignmentHole.solid = new G4IntersectionSolid(alignmentHole.name,
                                                  parent.solid,
                                                  holeCylinder,
                                                  colrot,
                                                  alignmentHole.centerInParent
                                                  );

    finishNesting(alignmentHole,
                  materialFinder.get("hall.insideMaterialName"),
                  0,
                  CLHEP::Hep3Vector(0,0,0),
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

    G4Polycone *plugCylinder = new G4Polycone(collimator.name()+"plugcomponent", 0, 2*M_PI,
                                              sizeof(zPlane)/sizeof(zPlane[0]),
                                              zPlane,
                                              rzero,
                                              rAlignmentPlug);

    VolumeInfo alignmentPlug(collimator.name()+"AlignmentPlug",
                             CLHEP::Hep3Vector(0,0,0),
                             alignmentHole.centerInWorld);

    alignmentPlug.solid = new G4IntersectionSolid(alignmentPlug.name,
                                                  parent.solid,
                                                  plugCylinder,
                                                  colrot,
                                                  alignmentHole.centerInParent
                                                  );

    finishNesting(alignmentPlug,
                  materialFinder.get("protonBeamDump.material.shielding"),
                  0,
                  CLHEP::Hep3Vector(0,0,0),
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
    // Collimator channel

    double rChannel[] = {
      collimator.channelRadius()[1],
      collimator.channelRadius()[1],
      collimator.channelRadius()[0],
      collimator.channelRadius()[0]
    };

    G4Polycone *channelCylinder = new G4Polycone(collimator.name()+"channel", 0, 2*M_PI,
                                                 sizeof(zPlane)/sizeof(zPlane[0]),
                                                 zPlane,
                                                 rzero,
                                                 rChannel);

    VolumeInfo channel(collimator.name()+"Channel",
                       CLHEP::Hep3Vector(0,0,0),
                       alignmentHole.centerInWorld);

    channel.solid = new G4IntersectionSolid(channel.name,
                                            parent.solid,
                                            channelCylinder,
                                            colrot,
                                            alignmentHole.centerInParent
                                            );

    finishNesting(channel,
                  materialFinder.get("hall.insideMaterialName"),
                  0,
                  CLHEP::Hep3Vector(0,0,0),
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
  void constructExtMonFNALMagnet(const ExtMonFNALMagnet& mag,
                                 const VolumeInfo& parent,
                                 const std::string& volNameSuffix,
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
      reg.add(mag.magnetRotationInMu2e().inverse()*parentRotationInMu2e);

    const CLHEP::Hep3Vector nx(1, 0, 0);
    const CLHEP::Hep3Vector ny(0, 1, 0);
    const CLHEP::Hep3Vector nz(0, 0, 1);
    AGDEBUG("magnetRotationInParent.inv * Nx = \n"
            <<*magnetRotationInParentInv*nx<<"\n"
            <<*magnetRotationInParentInv*ny<<"\n"
            <<*magnetRotationInParentInv*nz<<"\n"
            );

    AGDEBUG("mag.refPointInMu2e() = "<<mag.refPointInMu2e());
    AGDEBUG("mag.geometricCenterInMu2e() = "<<mag.geometricCenterInMu2e());
    AGDEBUG("magnet parent.centerInMu2e() = "<<parent.centerInMu2e()<<", in world = "<<parent.centerInWorld<<", in parent = "<<parent.centerInParent);
    AGDEBUG("magnet center in parent = "<<parentRotationInMu2e.inverse()*(mag.geometricCenterInMu2e() - parent.centerInMu2e()));

    const VolumeInfo magnetIron = nestBox("ExtMonFNAL"+volNameSuffix+"MagnetIron",
                                          mag.outerHalfSize(),
                                          materialFinder.get("extMonFNAL."+volNameSuffix+".magnet.material"),
                                          magnetRotationInParentInv,
                                          parentRotationInMu2e.inverse()*(mag.geometricCenterInMu2e() - parent.centerInMu2e()),
                                          parent, 0,
                                          config.getBool("extMonFNAL."+volNameSuffix+".magnet.iron.visible"),
                                          G4Colour::Magenta(),
                                          config.getBool("extMonFNAL."+volNameSuffix+".magnet.iron.solid"),
                                          forceAuxEdgeVisible,
                                          placePV,
                                          doSurfaceCheck
                                          );

    std::vector<double> apertureHalfSize(3);
    apertureHalfSize[0] = 0.5*mag.apertureWidth();
    apertureHalfSize[1] = 0.5*mag.apertureHeight();
    apertureHalfSize[2] = mag.magneticLength()/2;

    VolumeInfo magnetAperture = nestBox("ExtMonFNAL"+volNameSuffix+"MagnetAperture",
                                        apertureHalfSize,
                                        materialFinder.get("hall.insideMaterialName"),
                                        0,
                                        CLHEP::Hep3Vector(0, 0, 0),
                                        magnetIron.logical, 0,
                                        config.getBool("extMonFNAL."+volNameSuffix+".magnet.aperture.visible"),
                                        G4Colour::Grey(),
                                        config.getBool("extMonFNAL."+volNameSuffix+".magnet.aperture.solid"),
                                        forceAuxEdgeVisible,
                                        placePV,
                                        doSurfaceCheck
                                        );

    // The non-magnetic "margins" of the aperture to account for the difference
    // of physical and magnetic lengths
    if(mag.outerHalfSize()[2] - mag.magneticLength()/2 < 0) {
      throw cet::exception("GEOM")
        << "ExtMon magnet length/2 = "<<mag.outerHalfSize()[2]<<" < magnetic length/2 = "<< mag.magneticLength()/2 <<"for "<<volNameSuffix<<"\n";
    }

    std::vector<double> apertureMarginHalfSize(3);
    apertureMarginHalfSize[0] = 0.5*mag.apertureWidth();
    apertureMarginHalfSize[1] = 0.5*mag.apertureHeight();
    apertureMarginHalfSize[2] = (mag.outerHalfSize()[2] - mag.magneticLength()/2)/2;

    const double apertureMarginOffset = (mag.magneticLength()/2 + mag.outerHalfSize()[2])/2;

    VolumeInfo apertureMarginUp =
      nestBox("ExtMonFNAL"+volNameSuffix+"MagnetApertureMarginUp",
              apertureMarginHalfSize,
              materialFinder.get("hall.insideMaterialName"),
              0,
              CLHEP::Hep3Vector(0, 0, +apertureMarginOffset),
              magnetIron.logical, 0,
              config.getBool("extMonFNAL."+volNameSuffix+".magnet.aperture.visible"),
              G4Colour::Grey(),
              config.getBool("extMonFNAL."+volNameSuffix+".magnet.aperture.solid"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );

    VolumeInfo apertureMarginDn = 
      nestBox("ExtMonFNAL"+volNameSuffix+"MagnetApertureMarginDn",
              apertureMarginHalfSize,
              materialFinder.get("hall.insideMaterialName"),
              0,
              CLHEP::Hep3Vector(0, 0, -apertureMarginOffset),
              magnetIron.logical, 0,
              config.getBool("extMonFNAL."+volNameSuffix+".magnet.aperture.visible"),
              G4Colour::Grey(),
              config.getBool("extMonFNAL."+volNameSuffix+".magnet.aperture.solid"),
              forceAuxEdgeVisible,
              placePV,
              doSurfaceCheck
              );


    //----------------------------------------------------------------
    // Define the field in the magnet

    AGDEBUG("ExtMonFNAL "+volNameSuffix+" magnet field = "<<mag.bfield());

    G4MagneticField *field = reg.add(new G4UniformMagField(mag.bfield()));

    G4Mag_UsualEqRhs *rhs  = reg.add(new G4Mag_UsualEqRhs(field));

    G4MagIntegratorStepper *integrator = reg.add(new G4ExactHelixStepper(rhs));
    //G4MagIntegratorStepper *integrator = reg.add(new G4NystromRK4(rhs));

    const double stepMinimum = config.getDouble("extMonFNAL."+volNameSuffix+".magnet.stepMinimum", 1.0e-2 * CLHEP::mm /*The default from G4ChordFinder.hh*/);
    G4ChordFinder          *chordFinder = reg.add(new G4ChordFinder(field, stepMinimum, integrator));

    const double deltaOld = chordFinder->GetDeltaChord();
    chordFinder->SetDeltaChord(config.getDouble("extMonFNAL."+volNameSuffix+".magnet.deltaChord", deltaOld));
    AGDEBUG("chordFinder: using deltaChord = "<<chordFinder->GetDeltaChord()<<" (default = "<<deltaOld<<")");

    G4FieldManager *manager = reg.add(new G4FieldManager(field, chordFinder));

    AGDEBUG("orig: manager epsMin = "<<manager->GetMinimumEpsilonStep()
            <<", epsMax = "<<manager->GetMaximumEpsilonStep()
            <<", deltaOneStep = "<<manager->GetDeltaOneStep()
            );

    manager->SetMinimumEpsilonStep(config.getDouble("extMonFNAL."+volNameSuffix+".magnet.minEpsilonStep", manager->GetMinimumEpsilonStep()));
    manager->SetMaximumEpsilonStep(config.getDouble("extMonFNAL."+volNameSuffix+".magnet.maxEpsilonStep", manager->GetMaximumEpsilonStep()));
    manager->SetDeltaOneStep(config.getDouble("extMonFNAL."+volNameSuffix+".magnet.deltaOneStep", manager->GetDeltaOneStep()));

    AGDEBUG("new:  manager epsMin = "<<manager->GetMinimumEpsilonStep()
            <<", epsMax = "<<manager->GetMaximumEpsilonStep()
            <<", deltaOneStep = "<<manager->GetDeltaOneStep()
            );


    magnetIron.logical->SetFieldManager(manager, true);
    // No field in the margins
    apertureMarginUp.logical->SetFieldManager(0, true);
    apertureMarginDn.logical->SetFieldManager(0, true);
  }

  //================================================================
  void constructExtMonFNALBuilding(const VolumeInfo& collimator1Parent,
                                         const CLHEP::HepRotation& collimator1ParentRotationInMu2e,
                                         const VolumeInfo& mainParent,
                                         const CLHEP::HepRotation& mainParentRotationInMu2e,
                                         const SimpleConfig& config)
  {
    bool const forceAuxEdgeVisible = config.getBool("g4.forceAuxEdgeVisible",false);
    bool const placePV             = true;

    MaterialFinder materialFinder(config);

    GeomHandle<ProtonBeamDump> dump;
    GeomHandle<ExtMonFNALBuilding> emfb;

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();

    //----------------------------------------------------------------
    // Collimator2 shielding

    static CLHEP::HepRotation collimator2ParentRotationInMu2e(CLHEP::HepRotation::IDENTITY);
    collimator2ParentRotationInMu2e.rotateX( 90*CLHEP::degree);
    collimator2ParentRotationInMu2e.rotateZ( 90*CLHEP::degree);

    VolumeInfo coll2Shielding("ExtMonFNALColl2Shielding",
                            CLHEP::Hep3Vector(0, (emfb->roomInsideYmin()+emfb->roomInsideYmax())/2.0, 0)
                            - mainParent.centerInMu2e(),
                            mainParent.centerInWorld);

    coll2Shielding.solid = new G4ExtrudedSolid(coll2Shielding.name, emfb->coll2ShieldingOutline(),
                                             emfb->roomInsideFullHeight()/2.0,
                                             G4TwoVector(0,0), 1., G4TwoVector(0,0), 1.);

    finishNesting(coll2Shielding,
                  materialFinder.get("extMonFNAL.room.wall.materialName"),
                  &collimator2ParentRotationInMu2e,
                  coll2Shielding.centerInParent,
                  mainParent.logical,
                  0,
                  geomOptions->isVisible( "coll2Shielding" ),
                  G4Colour::Red() ,
                    geomOptions->isSolid( "coll2Shielding" ),
                    geomOptions->forceAuxEdgeVisible( "coll2Shielding" ),
                    geomOptions->placePV( "coll2Shielding" ),
                    geomOptions->doSurfaceCheck( "coll2Shielding" )
                  );

    //----------------------------------------------------------------
    // The filter channel

    AGDEBUG("emfb->collimator1CenterInMu2e() = "<<emfb->collimator1CenterInMu2e());
    AGDEBUG("collimator1Parent.centerInMu2e() = "<<collimator1Parent.centerInMu2e());

    constructCollimatorExtMonFNAL(emfb->collimator1(),
                                  collimator1Parent,
                                  collimator1ParentRotationInMu2e*(emfb->collimator1CenterInMu2e()-collimator1Parent.centerInMu2e()),
                                  collimator1ParentRotationInMu2e*emfb->collimator1RotationInMu2e(),
                                  config);

    constructExtMonFNALMagnet(emfb->filterMagnet(), mainParent, "filter", mainParentRotationInMu2e, config);

    constructCollimatorExtMonFNAL(emfb->collimator2(),
                                  coll2Shielding,
                                  collimator2ParentRotationInMu2e*(emfb->collimator2CenterInMu2e() - coll2Shielding.centerInMu2e()),
                                  collimator2ParentRotationInMu2e*emfb->collimator2RotationInMu2e(),
                                  config);

    // Test
    if(false) {
      G4Helper* _helper = &(*(art::ServiceHandle<G4Helper>()));
      const VolumeInfo& hall = _helper->locateVolInfo("HallAir");
      VolumeInfo test("emfMagnettest", emfb->filterMagnet().geometricCenterInMu2e() - hall.centerInMu2e(), hall.centerInWorld);
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

  } // constructExtMonFNALBuilding()
} // namespace mu2e
