// Andrei Gaponenko, 2012

#include "Mu2eG4/inc/constructExtMonFNAL.hh"

#include <iostream>
#include <iterator>
#include <algorithm>
#include <cmath>

#include "Geant4/G4Color.hh"
#include "Geant4/G4RotationMatrix.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4ExtrudedSolid.hh"
#include "Geant4/G4Trap.hh"
#include "Geant4/G4Orb.hh"
#include "Geant4/G4Box.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4ExtrudedSolid.hh"
#include "Geant4/G4IntersectionSolid.hh"
#include "Geant4/G4TwoVector.hh"
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

#include "GeometryService/inc/GeomHandle.hh"
#include "ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"

#include "Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Mu2eG4Helper/inc/AntiLeakRegistry.hh"
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
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNAL", "extMonFNAL" );
    geomOptions->loadEntry( config, "extMonFNAL"+collimator.name()+"alignmentHole", "extMonFNAL."+collimator.name()+".alignmentHole" );
    geomOptions->loadEntry( config, "extMonFNAL"+collimator.name()+"alignmentPlug", "extMonFNAL."+collimator.name()+".alignmentPlug" );
    geomOptions->loadEntry( config, "extMonFNAL"+collimator.name()+"channel",       "extMonFNAL."+collimator.name()+".channel" );

    const bool forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("extMonFNAL");
    const bool doSurfaceCheck       = geomOptions->doSurfaceCheck("extMonFNAL");
    const bool placePV              = geomOptions->placePV("extMonFNAL");


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
                  geomOptions->isVisible("extMonFNAL"+collimator.name()+"alignmentHole"),
                  G4Colour::Red(),//                  G4Colour::Cyan(),
                  geomOptions->isSolid("extMonFNAL"+collimator.name()+"alignmentHole"),
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
                  geomOptions->isVisible("extMonFNAL"+collimator.name()+"alignmentPlug"),
                  G4Colour(0.4, 0, 0), // G4Colour::Red(),
                  geomOptions->isSolid("extMonFNAL"+collimator.name()+"alignmentPlug"),
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
                  geomOptions->isVisible("extMonFNAL"+collimator.name()+"channel"),
                  G4Colour::Yellow(),
                  geomOptions->isSolid("extMonFNAL"+collimator.name()+"channel"),
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

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNAL", "extMonFNAL" );
    geomOptions->loadEntry( config, "extMonFNAL"+volNameSuffix+"magnetIron",     "extMonFNAL."+volNameSuffix+".magnet.iron" );
    geomOptions->loadEntry( config, "extMonFNAL"+volNameSuffix+"magnetAperture", "extMonFNAL."+volNameSuffix+".magnet.aperture" );

    const bool forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("extMonFNAL");
    const bool doSurfaceCheck       = geomOptions->doSurfaceCheck("extMonFNAL");
    const bool placePV              = geomOptions->placePV("extMonFNAL");

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

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
                                          geomOptions->isVisible("extMonFNAL"+volNameSuffix+"magnetIron"),
                                          G4Colour::Magenta(),
                                          geomOptions->isSolid("extMonFNAL"+volNameSuffix+"magnetIron"),
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
                                        geomOptions->isVisible("extMonFNAL"+volNameSuffix+"magnetAperture"),
                                        G4Colour::Grey(),
                                        geomOptions->isSolid("extMonFNAL"+volNameSuffix+"magnetAperture"),
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
              geomOptions->isVisible("extMonFNAL"+volNameSuffix+"magnetAperture"),
              G4Colour::Grey(),
              geomOptions->isSolid("extMonFNAL"+volNameSuffix+"magnetAperture"),
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
              geomOptions->isVisible("extMonFNAL"+volNameSuffix+"magnetAperture"),
              G4Colour::Grey(),
              geomOptions->isSolid("extMonFNAL"+volNameSuffix+"magnetAperture"),
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

    G4FieldManager *manager = new G4FieldManager(field, chordFinder);

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
    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNAL", "extMonFNAL" );
    geomOptions->loadEntry( config, "coll2Shielding",         "extMonFNAL.collimator2.shielding" );
    geomOptions->loadEntry( config, "coll2ShieldingHVACduct", "extMonFNAL.collimator2.shieldingHVACduct" );

    const bool forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("extMonFNAL");
    const bool doSurfaceCheck       = geomOptions->doSurfaceCheck("extMonFNAL");
    const bool placePV              = geomOptions->placePV("extMonFNAL");

    MaterialFinder materialFinder(config);

    GeomHandle<ProtonBeamDump> dump;
    GeomHandle<ExtMonFNALBuilding> emfb;


    static CLHEP::HepRotation shieldingRotationInMu2e = emfb->shieldingRotationInMu2e();
    const CLHEP::Hep3Vector offset(0.0,2*emfb->shieldingNHalfSize()[1],0.0);

    nestBox("ExtMonShieldingNorthsteel",
            emfb->shieldingNHalfSize(),
            materialFinder.get("protonBeamDump.material.core"),
            &shieldingRotationInMu2e,
            emfb->shieldingNCenterInMu2e() - mainParent.centerInMu2e(),
            mainParent.logical, 0,
            G4Colour::Red()
            );
    nestBox("ExtMonShieldingNorthconcrete",
            emfb->shieldingNHalfSize(),
            materialFinder.get("protonBeamDump.material.shielding"),
            &shieldingRotationInMu2e,
            emfb->shieldingNCenterInMu2e() - mainParent.centerInMu2e() + offset,
            mainParent.logical, 0,
            G4Colour::Grey()
            );
    nestBox("ExtMonShieldingSouthsteel",
            emfb->shieldingSHalfSize(),
            materialFinder.get("protonBeamDump.material.core"),
            &shieldingRotationInMu2e,
            emfb->shieldingSCenterInMu2e() - mainParent.centerInMu2e(),
            mainParent.logical, 0,
            G4Colour::Red()
            );
    nestBox("ExtMonShieldingSouthconcrete",
            emfb->shieldingSHalfSize(),
            materialFinder.get("protonBeamDump.material.shielding"),
            &shieldingRotationInMu2e,
            emfb->shieldingSCenterInMu2e() - mainParent.centerInMu2e() + offset,
            mainParent.logical, 0,
            G4Colour::Grey()
            );
    nestBox("ExtMonShieldingFloorConcrete",
            emfb->shieldingBHalfSize(),
            materialFinder.get("protonBeamDump.material.shielding"),
            &shieldingRotationInMu2e,
            emfb->shieldingBCenterInMu2e() - mainParent.centerInMu2e(),
            mainParent.logical, 0,
            G4Colour::Grey()
            );

    //----------------------------------------------------------------
    // Collimator2 shielding

    static CLHEP::HepRotation collimator2ParentRotationInMu2e = emfb->coll2ShieldingRotationInMu2e();

    VolumeInfo coll2Shielding("ExtMonFNALColl2Shielding",
                              emfb->coll2ShieldingCenterInMu2e() - mainParent.centerInMu2e(),
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
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );


    static const CLHEP::HepRotation HVACductRotInParent( collimator2ParentRotationInMu2e*emfb->shieldingRotationInMu2e().inverse() );

    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
    CLHEP::HepRotation *ductrot = reg.add(HVACductRotInParent.inverse());

    G4Tubs *holeCylinder = new G4Tubs( "holeCylinder", 0.0, emfb->HVACductRadius(), emfb->HVACductHalfLength(), 0.0, CLHEP::twopi );

    VolumeInfo HVACduct("coll2ShieldingHVACduct",
                        collimator2ParentRotationInMu2e*(emfb->HVACductCenterInMu2e() - emfb->coll2ShieldingCenterInMu2e()),
                        coll2Shielding.centerInWorld);

    HVACduct.solid = new G4IntersectionSolid(HVACduct.name,
                                             coll2Shielding.solid,
                                             holeCylinder,
                                             ductrot,
                                             HVACduct.centerInParent
                                             );

    finishNesting(HVACduct,
                  materialFinder.get("hall.insideMaterialName"),
                  0,
                  CLHEP::Hep3Vector(0,0,0),
                  coll2Shielding.logical,
                  0,
                  geomOptions->isVisible( "coll2ShieldingHVACduct" ),
                  G4Colour::G4Colour::Cyan(),
                  geomOptions->isSolid( "coll2ShieldingHVACduct" ),
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
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
    if (false) {
      Mu2eG4Helper* _helper = &(*(art::ServiceHandle<Mu2eG4Helper>()));
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
