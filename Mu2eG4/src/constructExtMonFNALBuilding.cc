// Andrei Gaponenko, 2012

#include "Offline/Mu2eG4/inc/constructExtMonFNAL.hh"

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
#include "Geant4/G4SubtractionSolid.hh"
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

#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ProtonBeamDumpGeom/inc/ProtonBeamDump.hh"
#include "Offline/ExtinctionMonitorFNAL/Geometry/inc/ExtMonFNALBuilding.hh"

#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/Mu2eG4Helper/inc/AntiLeakRegistry.hh"
#include "Offline/ConfigTools/inc/SimpleConfig.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/MaterialFinder.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/SensitiveDetectorName.hh"
#include "Offline/Mu2eG4/inc/FieldMgr.hh"
#include "Offline/GeomPrimitives/inc/TubsParams.hh"
#include "Offline/GeomPrimitives/inc/PolyconsParams.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/nestPolycone.hh"

// FIXME: should not need that
#include "Offline/GeometryService/inc/WorldG4.hh"

//#define AGDEBUG(stuff) std::cerr<<"AG: "<<__FILE__<<", line "<<__LINE__<<": "<<stuff<<std::endl;
#define AGDEBUG(stuff)

using namespace std;

namespace {

  double setDiameter(const vector<CLHEP::Hep2Vector>& v) {
    double minX = 0;
    double minY = 0;
    double maxX = 0;
    double maxY = 0;

    for (const auto& elem : v) {
        double currentX = elem.x();
        double currentY = elem.y();
        minX = (currentX <= minX) ? currentX : minX;
        minY = (currentY <= minY) ? currentY : minY;
        maxX = (currentX >= maxX) ? currentX : maxX;
        maxY = (currentY >= maxY) ? currentY : maxY;
    }

    return sqrt(pow((maxX-minX),2) + pow((maxY-minY),2));
  }
}

namespace mu2e {

  //================================================================
  void constructCollimatorExtMonFNAL(const ExtMonFNALBuilding::CollimatorExtMonFNAL& collimator,
                                     const VolumeInfo& parent,
                                     const CLHEP::Hep3Vector& collimatorCenterInParent,
                                     const CLHEP::HepRotation& collimatorRotationInParent,
                                     const SimpleConfig& config
                                     )
  {
    MaterialFinder materialFinder(config);
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();

    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( config, "extMonFNAL", "extMonFNAL" );
    geomOptions->loadEntry( config, "extMonFNAL"+collimator.name()+"channel", "extMonFNAL."+collimator.name()+".channel" );

    const bool forceAuxEdgeVisible  = geomOptions->forceAuxEdgeVisible("extMonFNAL");
    const bool doSurfaceCheck       = geomOptions->doSurfaceCheck("extMonFNAL");
    const bool placePV              = geomOptions->placePV("extMonFNAL");

    // The G4 interface we use through finishNesting() applies
    // the backwards interpretation of rotations.  Be consistent
    // and use the "backwards" interface of boolean solids.
    CLHEP::HepRotation *colrot = reg.add(collimatorRotationInParent.inverse());
    G4Material*  airMaterial = materialFinder.get("hall.insideMaterialName");


    //--------------------------------------------------------------------
    // Collimator Mother

    TubsParams motherParams{ 0, collimator.shotLinerOuterRadius(), 0.5*collimator.length() };

    VolumeInfo collimatorMother = nestTubs(collimator.name()+"Mother",
                                           motherParams,
                                           airMaterial,
                                           colrot,
                                           collimatorCenterInParent,
                                           parent,
                                           0,
                                           geomOptions->isVisible(collimator.name()+"Mother"),
                                           G4Colour::Red(),
                                           geomOptions->isSolid(collimator.name()+"Mother"),
                                           forceAuxEdgeVisible,
                                           placePV,
                                           doSurfaceCheck
                                           );

    //--------------------------------------------------------------------
    //parameters for nestPolycone

    // Fixme: protect against two planes at the same value of z.
    double const epsilon= ( collimator.radiusTransitiondZ() == 0. ) ? 0.1 : 0.;

    vector<double>zPlanes = {
      -0.5*collimator.length(),
      -0.5*collimator.radiusTransitiondZ()-epsilon,
      +0.5*collimator.radiusTransitiondZ()+epsilon,
      +0.5*collimator.length()
    };

    vector<double>rInner (zPlanes.size(), 0. );

    //--------------------------------------------------------------------
    // Outer Shot Liner

    TubsParams shotLinerOuterParams(collimator.shotLinerOuterRadius()
                                    -collimator.shotLinerOuterThickness(),
                                    collimator.shotLinerOuterRadius(),
                                    0.5*collimator.length()
                                    );

    VolumeInfo shotLinerOuter = nestTubs(collimator.name()+"shotLinerOuter",
                                         shotLinerOuterParams,
                                         findMaterialOrThrow("MildSteel"),
                                         0,
                                         CLHEP::Hep3Vector(0,0,0),
                                         collimatorMother,
                                         0,
                                         geomOptions->isVisible("extMonFNAL"+collimator.name()+"shotLinerOuter"),
                                         G4Colour::Red(),
                                         geomOptions->isSolid("extMonFNAL"+collimator.name()+"shotLinerOuter"),
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    //--------------------------------------------------------------------
    // Steel Shot

    vector<double> rInnerSteelShot = {
      collimator.shotLinerInnerRadius()[1] + collimator.shotLinerInnerThickness()[1],
      collimator.shotLinerInnerRadius()[1] + collimator.shotLinerInnerThickness()[1],
      collimator.shotLinerInnerRadius()[0] + collimator.shotLinerInnerThickness()[0],
      collimator.shotLinerInnerRadius()[0] + collimator.shotLinerInnerThickness()[0]
    };

    vector<double>rOuterSteelShot (zPlanes.size(),
                                   collimator.shotLinerOuterRadius()
                                   -collimator.shotLinerOuterThickness()
                                   );

    assert( zPlanes.size() == rInnerSteelShot.size() );

    PolyconsParams steelShotParams ( zPlanes, rInnerSteelShot, rOuterSteelShot );

    VolumeInfo steelShot = nestPolycone(collimator.name()+"SteelShot",
                                            steelShotParams,
                                            materialFinder.get("extMonFNAL."+collimator.name()+".shot.materialName"),
                                            0,
                                            CLHEP::Hep3Vector(0,0,0),
                                            collimatorMother,
                                            0,
                                            geomOptions->isVisible("extMonFNAL"+collimator.name()+"steelShot"),
                                            G4Colour::Red(),
                                            geomOptions->isSolid("extMonFNAL"+collimator.name()+"steelShot"),
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );

    //--------------------------------------------------------------------
    // Inner Shot Liner


    vector<double> rInnerInnerShotLinerRadius = {
      collimator.shotLinerInnerRadius()[1],
      collimator.shotLinerInnerRadius()[1],
      collimator.shotLinerInnerRadius()[0],
      collimator.shotLinerInnerRadius()[0]
    };

    vector<double> rOuterInnerShotLinerRadius = rInnerSteelShot;

    assert( zPlanes.size() == rOuterInnerShotLinerRadius.size() );

    PolyconsParams innerShotLinerParams ( zPlanes,
                                          rInnerInnerShotLinerRadius,
                                          rOuterInnerShotLinerRadius
                                          );

    VolumeInfo innerShotLiner = nestPolycone(collimator.name()+"InnerShotLiner",
                                             innerShotLinerParams,
                                             findMaterialOrThrow("MildSteel"),
                                             0,
                                             CLHEP::Hep3Vector(0,0,0),
                                             collimatorMother,
                                             0,
                                             geomOptions->isVisible("extMonFNAL"+collimator.name()+"innerShotLiner"),
                                             G4Colour::Red(),
                                             geomOptions->isSolid("extMonFNAL"+collimator.name()+"innerShotLiner"),
                                             forceAuxEdgeVisible,
                                             placePV,
                                             doSurfaceCheck
                                             );

    //--------------------------------------------------------------------
    // Alignment Plug Outer Shell

    vector<double> rOuterAlignmentPlugOuterShell = {
      collimator.alignmentPlugRadius()[1],
      collimator.alignmentPlugRadius()[1],
      collimator.alignmentPlugRadius()[0],
      collimator.alignmentPlugRadius()[0]
    };

    vector<double> rInnerAlignmentPlugOuterShell = {
      collimator.alignmentPlugRadius()[1] - collimator.alignmentPlugOuterShellThickness()[1],
      collimator.alignmentPlugRadius()[1] - collimator.alignmentPlugOuterShellThickness()[1],
      collimator.alignmentPlugRadius()[0] - collimator.alignmentPlugOuterShellThickness()[0],
      collimator.alignmentPlugRadius()[0] - collimator.alignmentPlugOuterShellThickness()[0]
    };

    assert( zPlanes.size() == rOuterAlignmentPlugOuterShell.size() );

    PolyconsParams alignmentPlugOuterShellParams ( zPlanes,
                                                   rInnerAlignmentPlugOuterShell,
                                                   rOuterAlignmentPlugOuterShell
                                                   );

    VolumeInfo alignmentPlugOuterShell = nestPolycone(collimator.name()+"AlignmentPlugOuterShell",
                                                      alignmentPlugOuterShellParams,
                                                      findMaterialOrThrow("MildSteel"),
                                                      0,
                                                      CLHEP::Hep3Vector(0,0,0),
                                                      collimatorMother,
                                                      0,
                                                      geomOptions->isVisible("extMonFNAL"+collimator.name()+"alignmentPlugOuterShell"),
                                                      G4Colour::Red(),
                                                      geomOptions->isSolid("extMonFNAL"+collimator.name()+"alignmentPlugOuterShell"),
                                                      forceAuxEdgeVisible,
                                                      placePV,
                                                      doSurfaceCheck
                                                      );

    //--------------------------------------------------------------------
    // Alignment Plug Concrete

    vector<double> rOuterAlignmentPlugConcrete = {
      collimator.alignmentPlugRadius()[1] - collimator.alignmentPlugOuterShellThickness()[1],
      collimator.alignmentPlugRadius()[1] - collimator.alignmentPlugOuterShellThickness()[1],
      collimator.alignmentPlugRadius()[0] - collimator.alignmentPlugOuterShellThickness()[0],
      collimator.alignmentPlugRadius()[0] - collimator.alignmentPlugOuterShellThickness()[0]
    };


    vector<double> rInnerAlignmentPlugConcrete = {
      collimator.channelRadius()[1] + collimator.alignmentPlugInnerShellThickness()[1],
      collimator.channelRadius()[1] + collimator.alignmentPlugInnerShellThickness()[1],
      collimator.channelRadius()[0] + collimator.alignmentPlugInnerShellThickness()[0],
      collimator.channelRadius()[0] + collimator.alignmentPlugInnerShellThickness()[0]
    };

    assert( zPlanes.size() == rOuterAlignmentPlugConcrete.size() );

    PolyconsParams alignmentPlugConcreteParams ( zPlanes,
                                                 rInnerAlignmentPlugConcrete,
                                                 rOuterAlignmentPlugConcrete
                                                 );

    VolumeInfo alignmentPlugConcrete = nestPolycone(collimator.name()+"AlignmentPlugConcrete",
                                                    alignmentPlugConcreteParams,
                                                    materialFinder.get("protonBeamDump.material.shielding"),
                                                    0,
                                                    CLHEP::Hep3Vector(0,0,0),
                                                    collimatorMother,
                                                    0,
                                                    geomOptions->isVisible("extMonFNAL"+collimator.name()+"alignmentPlugConcrete"),
                                                    G4Colour::Red(),
                                                    geomOptions->isSolid("extMonFNAL"+collimator.name()+"alignmentPlugConcrete"),
                                                    forceAuxEdgeVisible,
                                                    placePV,
                                                    doSurfaceCheck
                                                    );
    //--------------------------------------------------------------------
    // Alignment Plug Inner Shell

    vector<double> rOuterAlignmentPlugInnerShell = {
      collimator.channelRadius()[1] + collimator.alignmentPlugInnerShellThickness()[1],
      collimator.channelRadius()[1] + collimator.alignmentPlugInnerShellThickness()[1],
      collimator.channelRadius()[0] + collimator.alignmentPlugInnerShellThickness()[0],
      collimator.channelRadius()[0] + collimator.alignmentPlugInnerShellThickness()[0]
    };

    vector<double> rInnerAlignmentPlugInnerShell = {
      collimator.channelRadius()[1],
      collimator.channelRadius()[1],
      collimator.channelRadius()[0],
      collimator.channelRadius()[0]
    };

    assert( zPlanes.size() == rOuterAlignmentPlugInnerShell.size() );

    PolyconsParams alignmentPlugInnerShellParams ( zPlanes,
                                                   rInnerAlignmentPlugInnerShell,
                                                   rOuterAlignmentPlugInnerShell
                                                   );

    VolumeInfo alignmentPlugInnerShell = nestPolycone(collimator.name()+"AlignmentPlugInnerShell",
                                                      alignmentPlugInnerShellParams,
                                                      findMaterialOrThrow("MildSteel"),
                                                      0,
                                                      CLHEP::Hep3Vector(0,0,0),
                                                      collimatorMother,
                                                      0,
                                                      geomOptions->isVisible("extMonFNAL"+collimator.name()+"alignmentPlugInnerShell"),
                                                      G4Colour::Red(),
                                                      geomOptions->isSolid("extMonFNAL"+collimator.name()+"alignmentPlugInnerShell"),
                                                      forceAuxEdgeVisible,
                                                      placePV,
                                                      doSurfaceCheck
                                                      );


    //--------------------------------------------------------------------
    // Collimator channel

    vector<double> rOuterChannel = {
      collimator.channelRadius()[1],
      collimator.channelRadius()[1],
      collimator.channelRadius()[0],
      collimator.channelRadius()[0]
    };

    assert( zPlanes.size() == rOuterChannel.size() );

    PolyconsParams channelParams ( zPlanes, rInner, rOuterChannel );

    VolumeInfo channel = nestPolycone(collimator.name()+"Channel",
                                      channelParams,
                                      airMaterial,
                                      0,
                                      CLHEP::Hep3Vector(0,0,0),
                                      collimatorMother,
                                      0,
                                      geomOptions->isVisible("extMonFNAL"+collimator.name()+"channel"),
                                      G4Colour::Red(),
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
    CLHEP::HepRotation *magnetRotationInParentInv = reg.add(mag.magnetRotationInMu2e().inverse()*parentRotationInMu2e);

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
    AntiLeakRegistry& reg = art::ServiceHandle<Mu2eG4Helper>()->antiLeakRegistry();
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

    //--------------------------------------------------------------------
    // Collimator2 shielding

    CLHEP::HepRotation *pshieldingRot = reg.add(new CLHEP::HepRotation());
    CLHEP::HepRotation& shieldingRot = *pshieldingRot;
    pshieldingRot->rotateX( 90*CLHEP::degree);
    pshieldingRot->rotateZ( 90*CLHEP::degree);

    double subCylinderLength = setDiameter(emfb->coll2ShieldingOutline());

    CLHEP::Hep3Vector subCylOffsetInParent = shieldingRot *(emfb->collimator2CenterInMu2e() - emfb->coll2ShieldingCenterInMu2e());

    static CLHEP::HepRotation collimator2ParentRotationInMu2e = emfb->coll2ShieldingRotationInMu2e();

    CLHEP::HepRotation *subCylinderRotation = reg.add(emfb->collimator2RotationInMu2e().inverse()*collimator2ParentRotationInMu2e.inverse());

    G4Tubs* subCylinder = new G4Tubs("ExtMonFNALCollimator2Hole",
                                     0.*CLHEP::mm,
                                     emfb->collimator2().shotLinerOuterRadius(),
                                     subCylinderLength,
                                     0,
                                     CLHEP::twopi
                                     );

    G4ExtrudedSolid* coll2ShieldingExtrusion= new G4ExtrudedSolid( "ExtMonFNALCollimator2ShieldingExtrusion",
                                                                   emfb->coll2ShieldingOutline(),
                                                                   emfb->roomInsideFullHeight()/2.0,
                                                                   G4TwoVector(0,0),
                                                                   1.,
                                                                   G4TwoVector(0,0), 1.
                                                                 );

    VolumeInfo coll2Shielding("ExtMonFNALColl2Shielding",
                              emfb->coll2ShieldingCenterInMu2e() - mainParent.centerInMu2e(),
                              mainParent.centerInWorld
                              );

    coll2Shielding.solid = new G4SubtractionSolid(coll2Shielding.name,
                                                  coll2ShieldingExtrusion,
                                                  subCylinder,
                                                  subCylinderRotation,
                                                  subCylOffsetInParent
                                                  );

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

    //--------------------------------------------------------------------
    // The filter channel

    AGDEBUG("emfb->collimator1CenterInMu2e() = "<<emfb->collimator1CenterInMu2e());
    AGDEBUG("collimator1Parent.centerInMu2e() = "<<collimatomaterialFinder.getr1Parent.centerInMu2e());


    constructCollimatorExtMonFNAL(emfb->collimator1(),
                                  mainParent,
                                  emfb->collimator1CenterInMu2e()-mainParent.centerInMu2e(),
                                  emfb->collimator1RotationInMu2e(),
                                  config);

    constructExtMonFNALMagnet(emfb->filterMagnet(), mainParent, "filter", mainParentRotationInMu2e, config);

    constructCollimatorExtMonFNAL(emfb->collimator2(),
                                  mainParent,
                                  emfb->collimator2CenterInMu2e() - mainParent.centerInMu2e(),
                                  emfb->collimator2RotationInMu2e(),
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

    //--------------------------------------------------------------------

  } // constructExtMonFNALBuilding()
} // namespace mu2e
