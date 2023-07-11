//
// Free function to create  Production Solenoid and Production Target.
//
//
// Original author KLG based on Mu2eWorld constructPS
//
// Notes:
// Construct the PS. Parent volume is the air inside of the hall.

// C++ includes
#include <iostream>

// Mu2e includes.
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/ProductionSolenoidGeom/inc/ProductionSolenoid.hh"
#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/GeomPrimitives/inc/Polycone.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/ProductionTargetGeom/inc/ProductionTarget.hh"
#include "Offline/ProductionSolenoidGeom/inc/PSEnclosure.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/constructPS.hh"
#include "Offline/Mu2eG4/inc/constructPSShield.hh"
#include "Offline/Mu2eG4/inc/constructTargetPS.hh"
#include "Offline/Mu2eG4/inc/constructHaymanRings.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"

#include "Offline/ProductionSolenoidGeom/inc/PSVacuum.hh"

// FIXME: only needed when the constructPSShield() call below is conditional
#include "Offline/ProductionSolenoidGeom/inc/PSShield.hh"

// G4 includes
#include "Geant4/G4ThreeVector.hh"
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4SDManager.hh"
#include "Geant4/G4LogicalVolume.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4SubtractionSolid.hh"

#include "cetlib_except/exception.h"
#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  G4LogicalVolume* constructPS(VolumeInfo const & parent, SimpleConfig const & _config) {

    ProductionSolenoid const & psgh = *(GeomHandle<ProductionSolenoid>());

    Tube const & psVacVesselInnerParams     = *psgh.getVacVesselInnerParamsPtr();
    Tube const & psVacVesselOuterParams     = *psgh.getVacVesselOuterParamsPtr();

    // Extract some information from the config file.

    int verbosityLevel                  = _config.getInt("PS.verbosityLevel");

    std::string targetPS_model = _config.getString("targetPS_model");

    G4Material* psVacVesselMaterial = findMaterialOrThrow(psVacVesselInnerParams.materialName());

    verbosityLevel >0 &&
      cout << __func__ << " verbosityLevel                   : " << verbosityLevel  << endl;

    G4GeometryOptions* geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "PS", "PS" );

    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

    // we need to replace this in future getting/moving the info from/to geom service

    // VacVessel

    // VacVessel is a "G4Tubs"; or a combination of those

    // Build the barrel of the cryostat, the length is the "outer length"

    VolumeInfo psVacVesselInnerInfo =
      nestTubs("PSVacVesselInner",
               psVacVesselInnerParams.getTubsParams(),
               psVacVesselMaterial,
               0,
               psVacVesselInnerParams.originInMu2e()-_hallOriginInMu2e,
               parent,
               0,
               G4Colour::Green(),
               "PS"
               );

    VolumeInfo psVacVesselOuterInfo =
      nestTubs("PSVacVesselOuter",
               psVacVesselOuterParams.getTubsParams(),
               psVacVesselMaterial,
               0,
               psVacVesselOuterParams.originInMu2e()-_hallOriginInMu2e,
               parent,
               0,
               G4Colour::Green(),
               "PS"
               );


    // Rings
    Tube const & psRing1Params = *psgh.getRing1ParamsPtr();
    Tube const & psRing2Params = *psgh.getRing2ParamsPtr();
    G4Material* ringMaterial = findMaterialOrThrow(psRing1Params.materialName());
    VolumeInfo psRing1Info =
      nestTubs("PSRing1", psRing1Params.getTubsParams(),
               ringMaterial, 0,
               psRing1Params.originInMu2e() - _hallOriginInMu2e,
               parent, 0, G4Colour::Blue(), "PS");
    VolumeInfo psRing2Info =
      nestTubs("PSRing2", psRing2Params.getTubsParams(),
               ringMaterial, 0,
               psRing2Params.originInMu2e() - _hallOriginInMu2e,
               parent, 0, G4Colour::Blue(), "PS");

    // two endplates

    Tube const & psVacVesselEndPlateDParams = *psgh.getVacVesselEndPlateDParamsPtr();
    Tube const & psVacVesselEndPlateUParams = *psgh.getVacVesselEndPlateUParamsPtr();

    VolumeInfo psVacVesselEndPlateDInfo =
      nestTubs("PSVacVesselEndPlateD",
               psVacVesselEndPlateDParams.getTubsParams(),
               psVacVesselMaterial,
               0,
               psVacVesselEndPlateDParams.originInMu2e()-_hallOriginInMu2e,
               parent,
               0,
               G4Colour::Yellow(),
               "PS"
               );

    VolumeInfo psVacVesselEndPlateUInfo
      = nestTubs("PSVacVesselEndPlateU",
                 psVacVesselEndPlateUParams.getTubsParams(),
                 psVacVesselMaterial,
                 0,
                 psVacVesselEndPlateUParams.originInMu2e()-_hallOriginInMu2e,
                 parent,
                 0,
                 G4Colour::Red(),
                 "PS"
                 );



    // Put vacuum inside vacuum vessel
    double vacRIn = psVacVesselInnerParams.getTubsParams().outerRadius();
    double vacROut = psVacVesselOuterParams.getTubsParams().innerRadius();
    double vacHalfLength = psVacVesselOuterParams.getTubsParams().zHalfLength() - 2.*psVacVesselEndPlateUParams.getTubsParams().zHalfLength();

    VolumeInfo psVacuumVesselVacuumInfo
      = nestTubs ( "psVacuumVesselVacuum",
                   TubsParams(vacRIn, vacROut, vacHalfLength),
                   findMaterialOrThrow("PSVacuum"),
                   0,
                   psVacVesselOuterParams.originInMu2e() - _hallOriginInMu2e,
                   parent, 0, G4Colour::White(),
                   "PS" );



    Polycone const & psCoilShellParams = *psgh.getCoilShellParamsPtr();

    //Coil "Outer Shell"

    //we will place the shell inside the parent, which used to be the hall,
    // but is now the vacuum in the PS cryo.

    string const psCoilShellName = "PSCoilShell";

    VolumeInfo psCoilShellInfo(psCoilShellName,
                               psCoilShellParams.originInMu2e()-psVacuumVesselVacuumInfo.centerInMu2e(),
                               psVacuumVesselVacuumInfo.centerInWorld);

    psCoilShellInfo.solid  =  new G4Polycone( psCoilShellName,
                                              psCoilShellParams.phi0(),
                                              psCoilShellParams.phiTotal(),
                                              psCoilShellParams.numZPlanes(),
                                              &psCoilShellParams.zPlanes()[0],
                                              &psCoilShellParams.rInner()[0],
                                              &psCoilShellParams.rOuter()[0]);

    G4Material* psCoilShellMaterial = findMaterialOrThrow(psCoilShellParams.materialName());

    finishNesting(psCoilShellInfo,
                  psCoilShellMaterial,
                  0,
                  psCoilShellParams.originInMu2e()-psVacuumVesselVacuumInfo.centerInMu2e(),
                  psVacuumVesselVacuumInfo.logical,
                  0,
                  G4Colour::White(),
                  "PS"
                  );

    // the superconducting Coils

    //the coils are "G4Tubs" placed inside the Coil Shell

    Tube const & psCoil1Params        = *psgh.getCoil1ParamsPtr();

    // all coil materials are the same...
    G4Material* psCoilMaterial        = findMaterialOrThrow(psCoil1Params.materialName());

    VolumeInfo psCoil1Info = nestTubs("PSCoil1",
                                      psCoil1Params.getTubsParams(),
                                      psCoilMaterial,
                                      0,
                                      psCoil1Params.originInMu2e()- psCoilShellParams.originInMu2e(),
                                      psCoilShellInfo,
                                      0,
                                      G4Colour::Red(),
                                      "PS"
                                      );


    Tube const & psCoil2Params        = *psgh.getCoil2ParamsPtr();

    VolumeInfo psCoil2Info = nestTubs("PSCoil2",
                                      psCoil2Params.getTubsParams(),
                                      psCoilMaterial,
                                      0,
                                      psCoil2Params.originInMu2e()- psCoilShellParams.originInMu2e(),
                                      psCoilShellInfo,
                                      0,
                                      G4Colour::Red(),
                                      "PS"
                                      );


    Tube const & psCoil3Params        = *psgh.getCoil3ParamsPtr();

    VolumeInfo psCoil3Info = nestTubs("PSCoil3",
                                      psCoil3Params.getTubsParams(),
                                      psCoilMaterial,
                                      0,
                                      psCoil3Params.originInMu2e()- psCoilShellParams.originInMu2e(),
                                      psCoilShellInfo,
                                      0,
                                      G4Colour::Red(),
                                      "PS"
                                      );


    if ( verbosityLevel > 0) {

      cout << __func__ << " parent.centerInMu2e()                     : "
           << parent.centerInMu2e() << endl;

      cout << __func__ << " psVacVesselInnerParams.originInMu2e()     : "
           << psVacVesselInnerParams.originInMu2e() << endl;
      cout << __func__ << " psVacVesselOuterParams.originInMu2e()     : "
           << psVacVesselOuterParams.originInMu2e() << endl;

      cout << __func__ << " psVacVesselEndPlateDParams.originInMu2e() : "
           << psVacVesselEndPlateDParams.originInMu2e() << endl;
      cout << __func__ << " psVacVesselEndPlateUParams.originInMu2e() : "
           << psVacVesselEndPlateUParams.originInMu2e() << endl;

      cout << __func__ << " psCoilShellParams.originInMu2e()          : "
           << psCoilShellParams.originInMu2e() << endl;

      cout << __func__ << " psCoilShellInfo.centerInMu2e()            : "
           << psCoilShellInfo.centerInMu2e() <<  endl;

      cout << __func__ << " psCoil1Params.originInMu2e()              : "
           << psCoil1Params.originInMu2e() << endl;
      cout << __func__ << " psCoil1 local Offset                      : "
           << psCoil1Params.originInMu2e()- psCoilShellInfo.centerInMu2e() << endl;

      cout << __func__ << " psCoil2Params.originInMu2e()              : "
           << psCoil2Params.originInMu2e() << endl;
      cout << __func__ << " psCoil2 local Offset                      : "
           << psCoil2Params.originInMu2e()- psCoilShellInfo.centerInMu2e() << endl;

      cout << __func__ << " psCoil3Params.originInMu2e()              : "
           << psCoil3Params.originInMu2e() << endl;
      cout << __func__ << " psCoil3 local Offset                      : "
           << psCoil3Params.originInMu2e()- psCoilShellInfo.centerInMu2e() << endl;

    }

    // Build the main PS vacuum body.

    // we may need to adjust it to make sure it conforms to the newer
    // drawings where the collimators enter the PS area


    Tube const & psVacuumParams  = GeomHandle<PSVacuum>()->vacuum();
    GeomHandle<PSEnclosure> pse;
    VolumeInfo psVacuumInfo("PSVacuum", psVacuumParams.originInMu2e() - _hallOriginInMu2e, parent.centerInWorld);
    if(pse->version() < 3) {
      psVacuumInfo = nestTubs( "PSVacuum",
                               psVacuumParams.getTubsParams(),
                               findMaterialOrThrow(psVacuumParams.materialName()),
                               0,
                               psVacuumParams.originInMu2e()-_hallOriginInMu2e,
                               parent,
                               0,
                               G4Colour::Blue(),
                               "PS"
                               );
    } else { //remove any of the PS vacuum from the PS end cap region
      CLHEP::Hep3Vector extraOffset(0.,0.,pse->getExtraOffset());
      const TubsParams& psVacuumTubsParams = psVacuumParams.getTubsParams();
      G4Tubs* psVacuumTubs = new G4Tubs("PSVacuumTube",
                                        psVacuumTubsParams.innerRadius(),
                                        psVacuumTubsParams.outerRadius(),
                                        psVacuumTubsParams.zHalfLength(),
                                        0., CLHEP::twopi);
      auto polycone = pse->endPlatePolycone();
      std::vector<double> rInners, rOuters(polycone.rOuter()), zPlanes(polycone.zPlanes());
      if ( verbosityLevel > 1) {
        cout << __func__ << " printing PS enclosure planes to subtract from PS vacuum:\n";
      }
      for(unsigned iplane = 0; iplane <= polycone.numZPlanes(); ++iplane) {
        rInners.push_back(0.); //clear from center to outer radius
        if( verbosityLevel > 1) {
          if(iplane < polycone.numZPlanes()) {
            cout << " (r1, r2, z) = (" << rInners[iplane] << ", " << rOuters[iplane]
                 << ", " << zPlanes[iplane] << ")\n";
          }
        }
      }
      rOuters.insert(rOuters.begin(), rOuters.front());
      zPlanes.insert(zPlanes.begin(), zPlanes.front()-1.e3); //clear beyond the endcap as well
      G4Polycone* basecone = new G4Polycone("PSEnclosureEndplateBaseCone_FromPS",
                                            polycone.phi0(),
                                            polycone.phiTotal(),
                                            zPlanes.size(),
                                            &zPlanes[0],
                                            &rInners[0],
                                            &rOuters[0]
                                            );
      G4SubtractionSolid* psVacuumSolid = new G4SubtractionSolid("PSVacuum", psVacuumTubs, basecone,
                                                                 nullptr, //assume no relative rotation
                                                                 polycone.originInMu2e() - psVacuumParams.originInMu2e() + extraOffset);
      psVacuumInfo.solid = psVacuumSolid;
      finishNesting(psVacuumInfo,
                    findMaterialOrThrow(psVacuumParams.materialName()),
                    nullptr,
                    psVacuumParams.originInMu2e()-_hallOriginInMu2e,
                    parent.logical,
                    0,
                    G4Colour::Blue(),
                    "PS"
                    );
    }
    //    std::cout << "inside " << __func__ << psVacuumParams.originInMu2e() << " " << _hallOriginInMu2e << std::endl;

    // Build the production target.

    if ( targetPS_model == "MDC2018" ){
      verbosityLevel> 0 && std::cout << __func__ << "MDC 2018 target" << std::endl;
      constructTargetPS(psVacuumInfo, _config );
    } else if ( targetPS_model == "HaymanLowerDensity" ){
      verbosityLevel> 0 && std::cout << __func__ << "HaymanLowerDensity target" << std::endl;
      constructHaymanRings(psVacuumInfo, _config);
    } else if (targetPS_model == "Hayman_v_2_0"){
      verbosityLevel> 0 && std::cout << __func__ << "Hayman 2.0 target" << std::endl;
      constructTargetPS(psVacuumInfo, _config );
    } else{
      throw cet::exception("CONFIG")
        << "In constructPS.cc unrecognized production target model name: "
        << targetPS_model
        << "\n";
    }

   // FIXME: make unconditional
    if(art::ServiceHandle<GeometryService>()->hasElement<PSShield>()) {
      constructPSShield(psVacuumInfo, _config);
    }


      return psVacuumInfo.logical;

  } // end Mu2eWorld::constructPS
}
