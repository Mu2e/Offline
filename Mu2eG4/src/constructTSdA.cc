//
// Free function to create TSdA Neutron Absorber in G4
//
//
// Original author KLG
//
// Notes:
// Construct the Internal Neutron Absorbers in G4
// David Norvil Brown (the other one):  rename to TSdA for consistency
// with TDR, and update - May 2015.

// art includes
#include "art/Framework/Services/Registry/ServiceDefinitionMacros.h"

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Units/PhysicalConstants.h"

// Mu2e includes.

#include "Offline/Mu2eG4/inc/constructTSdA.hh"
#include "Offline/BeamlineGeom/inc/Beamline.hh"
#include "Offline/BeamlineGeom/inc/StraightSection.hh"
#include "Offline/Mu2eG4Helper/inc/VolumeInfo.hh"
#include "Offline/GeometryService/inc/GeomHandle.hh"
#include "Offline/GeometryService/inc/GeometryService.hh"
#include "Offline/GeometryService/inc/G4GeometryOptions.hh"
#include "Offline/Mu2eG4Helper/inc/Mu2eG4Helper.hh"
#include "Offline/DetectorSolenoidGeom/inc/DetectorSolenoid.hh"
#include "Offline/BeamlineGeom/inc/TSdA.hh"
#include "Offline/Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Offline/Mu2eG4/inc/nestBox.hh"
#include "Offline/Mu2eG4/inc/nestTubs.hh"
#include "Offline/Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "Geant4/G4Material.hh"
#include "Geant4/G4Color.hh"
#include "Geant4/G4VSolid.hh"
#include "Geant4/G4Tubs.hh"
#include "Geant4/G4Cons.hh"
#include "Geant4/G4Polycone.hh"
#include "Geant4/G4VPhysicalVolume.hh"
#include "Geant4/G4SubtractionSolid.hh"

// C++ includes
#include <string>

using namespace std;

namespace mu2e {

  void constructTSdA(SimpleConfig const & _config){

    int const verbosityLevel = _config.getInt("tsda.verbosityLevel",0);
    double tmpRin = _config.getDouble("tsda.rin",-1.0)*CLHEP::mm;


    const auto geomOptions = art::ServiceHandle<GeometryService>()->geomOptions();
    geomOptions->loadEntry( _config, "tsda", "tsda");

    bool const NAVisible           = geomOptions->isVisible("tsda");
    bool const NASolid             = geomOptions->isSolid("tsda");
    const bool doSurfaceCheck      = geomOptions->doSurfaceCheck("tsda");
    const bool forceAuxEdgeVisible = geomOptions->forceAuxEdgeVisible("tsda");
    const bool placePV             = geomOptions->placePV("tsda");

    // Access to the Mu2eG4HelperService.
    Mu2eG4Helper* _helper = &(*(art::ServiceHandle<Mu2eG4Helper>()));

    // now constructing the internal neutron absorber
    // it is placed inside DS2Vacuum & DS3Vacuum like the protonabs1 & 2

    // Fetch TSdA geometry
    GeomHandle<TSdA> atsd;

    int version = atsd->version();
    if ( tmpRin < 0. ) {  // check for negative value indicating to use the default value
      // TS5 outer radius
      GeomHandle<Beamline> beamg;
      const StraightSection * ts5out = beamg->getTS().getTSCryo<StraightSection>( TransportSolenoid::TSRegion::TS5,
                                                                                  TransportSolenoid::TSRadialPart::IN );
      tmpRin = ts5out->rIn();
    }
    if ( verbosityLevel > 0) {
      cout << __func__ << " TSdA rin                      : "<< tmpRin << endl;
      cout << __func__ << " TSdA rout                          : "<< atsd->r4() << endl;
    }

    // we need to calculate where the DS2Vacuum volume is;
    // the int. neutron absorber must fit inside it

    // Fetch DS geom object
    GeomHandle<DetectorSolenoid> ds;

    double const ds2FrontZ = ds->vac_zLocDs23Split() - 2.*ds->vac_halfLengthDs2();
    if ( verbosityLevel > 0) {
      cout << __func__ << " DS2Vacuum extent               : ["
           << ds2FrontZ << " , "
           << ds->vac_zLocDs23Split() << endl;
    }

    // certain combinations are illegal; we will assume that the
    // conical absorber is fully contained in ds2


    // Get DS2Vacuum & DS3Vacuum info
    VolumeInfo const & ds2VacuumInfo = _helper->locateVolInfo("DS2Vacuum");
    VolumeInfo const & ds3VacuumInfo = _helper->locateVolInfo("DS3Vacuum");

    CLHEP::Hep3Vector ATSDOffset = atsd->position() - ds2VacuumInfo.centerInMu2e();

    CLHEP::Hep3Vector ATSD4Offset(0.0,0.0,-ds->vac_halfLengthDs2() + atsd->halfLength4());

    if ( verbosityLevel > 0) {
      cout << __func__ << " DS2VacuumInfo.centerInMu2e()  : " << ds2VacuumInfo.centerInMu2e() << endl;
      cout << __func__ << " DS3VacuumInfo.centerInMu2e()  : " << ds3VacuumInfo.centerInMu2e() << endl;
      cout << __func__ << " TSdA Offset                   : "  << ATSDOffset  << endl;
      cout << __func__ << " TSdA Offset2                       : "  << ATSD4Offset  << endl;
      cout << __func__ << " TSdA center  in Mu2e          : "  << (ds2VacuumInfo.centerInMu2e() + ATSDOffset ).z()  << endl;
      cout << __func__ << " TSdA center2 in Mu2e          : "  << (ds2VacuumInfo.centerInMu2e() + ATSD4Offset).z()  << endl;
    }

    TubsParams ATSD4Params( tmpRin, atsd->r4(), atsd->halfLength4() );

    G4Material* ATSD4Material = findMaterialOrThrow( atsd->material4() );

    VolumeInfo info;
    info.name     = "TSdA4";
    info.solid    = new G4Tubs("TSdA4",
                               ATSD4Params.data()[0], ATSD4Params.data()[1],
                               ATSD4Params.data()[2], ATSD4Params.data()[3],
                               ATSD4Params.data()[4]);

    // Optionally add a cut-out to the disk
    const bool add_cutout = _config.getBool("tsda.cutout.build", false);
    if(add_cutout) {
      const double r_in       = _config.getDouble("tsda.cutout.rin" )*CLHEP::mm;
      const double r_out      = _config.getDouble("tsda.cutout.rout")*CLHEP::mm;
      const double phi_0      = _config.getDouble("tsda.cutout.phi0")*CLHEP::twopi/360.;
      const double dphi       = _config.getDouble("tsda.cutout.dphi")*CLHEP::twopi/360.;
      const double x_0        = _config.getDouble("tsda.cutout.x")*CLHEP::mm;
      const double y_0        = _config.getDouble("tsda.cutout.y")*CLHEP::mm;

      auto cutout = new G4Tubs("TSdA4Cutout"
                               , r_in
                               , r_out
                               , atsd->halfLength4() * 1.01 // ensure it's bigger
                               , phi_0
                               , dphi
                               );

      auto remainingDisk = new G4SubtractionSolid("TSdA4RemainingDisk",
                                                  info.solid,
                                                  cutout,
                                                  nullptr,
                                                  G4ThreeVector(x_0, y_0, 0)
                                                  );

      info.solid = remainingDisk;
    }

    finishNesting(info,
                  ATSD4Material,
                  0,
                  (version <= 1) ? ATSD4Offset : ATSDOffset, //in old versions, continue overriding position
                  ds2VacuumInfo.logical,
                  0,
                  NAVisible,
                  G4Colour::Cyan(),
                  NASolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );


    // Construct an additional object if requested (e.g. additional shielding)
    const bool build_extra = _config.getBool("tsda.extra.build", false);
    if(build_extra) {
      const double r_in       = _config.getDouble("tsda.extra.rin")*CLHEP::mm;
      const double r_out      = _config.getDouble("tsda.extra.rout")*CLHEP::mm;
      const double halflength = _config.getDouble("tsda.extra.halflength")*CLHEP::mm;
      const double z0         = _config.getDouble("tsda.extra.z0")*CLHEP::mm;
      G4Material* material    = findMaterialOrThrow(_config.getString("tsda.extra.material"));
      TubsParams params(r_in, r_out, halflength);
      CLHEP::Hep3Vector offset(0.,0.,z0 - ds2VacuumInfo.centerInMu2e().z());
      nestTubs("TSdExtra",
               params,
               material,
               0,
               offset,
               ds2VacuumInfo,
               0,
               NAVisible,
               G4Colour::Cyan(),
               NASolid,
               forceAuxEdgeVisible,
               placePV,
               doSurfaceCheck
               );

    }

    // Optionally construct a vector of tubes
    const bool build_tubes = _config.getBool("tsda.tubes.build", false);
    if(build_tubes) {
      const size_t n_tubes = _config.getInt("tsda.tubes.n");
      std::vector<double>      tube_rins       ; _config.getVectorDouble("tsda.tubes.rins"       , tube_rins       , n_tubes);
      std::vector<double>      tube_routs      ; _config.getVectorDouble("tsda.tubes.routs"      , tube_routs      , n_tubes);
      std::vector<double>      tube_halflengths; _config.getVectorDouble("tsda.tubes.halflengths", tube_halflengths, n_tubes);
      std::vector<double>      tube_z0s        ; _config.getVectorDouble("tsda.tubes.z0s"        , tube_z0s        , n_tubes);
      std::vector<double>      tube_dxs        ; _config.getVectorDouble("tsda.tubes.dxs"        , tube_dxs        , n_tubes);
      std::vector<double>      tube_dys        ; _config.getVectorDouble("tsda.tubes.dys"        , tube_dys        , n_tubes);
      std::vector<std::string> tube_materials  ; _config.getVectorString("tsda.tubes.materials"  , tube_materials  , n_tubes);
      if(tube_rins.size() != n_tubes || tube_routs.size() != n_tubes || tube_halflengths.size() != n_tubes ||
         tube_z0s.size() != n_tubes || tube_materials.size() != n_tubes || tube_dxs.size() != n_tubes || tube_dys.size() != n_tubes) {
        throw std::runtime_error("Size of tube parameter vectors must match n_tubes");
      }

      // Add each tube to the geometry
      for(size_t i = 0; i < n_tubes; ++i) {
        TubsParams params(tube_rins[i]*CLHEP::mm, tube_routs[i]*CLHEP::mm, tube_halflengths[i]*CLHEP::mm);
        CLHEP::Hep3Vector offset(tube_dxs[i]*CLHEP::mm, tube_dys[i]*CLHEP::mm,
                                 tube_z0s[i]*CLHEP::mm - ds2VacuumInfo.centerInMu2e().z());
        G4Material* material = findMaterialOrThrow(tube_materials[i]);
        nestTubs(std::string("TSdATube") + std::to_string(i),
                 params,
                 material,
                 0,
                 offset,
                 ds2VacuumInfo,
                 0,
                 NAVisible,
                 G4Colour::Cyan(),
                 NASolid,
                 forceAuxEdgeVisible,
                 placePV,
                 doSurfaceCheck
                 );
      }
    } // end of optional tubes

  } // end of constructTSdA;

}
