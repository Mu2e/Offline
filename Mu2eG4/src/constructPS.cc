//
// Free function to create  Production Solenoid and Production Target.
//
// $Id: constructPS.cc,v 1.10 2012/03/10 00:00:26 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/10 00:00:26 $
//
// Original author KLG based on Mu2eWorld constructPS
//
// Notes:
// Construct the PS. Parent volume is the air inside of the hall.

// C++ includes
#include <iostream>

// Mu2e includes.
#include "BeamlineGeom/inc/Beamline.hh"
#include "G4Helper/inc/VolumeInfo.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "ProductionTargetGeom/inc/ProductionTarget.hh"
#include "Mu2eG4/inc/findMaterialOrThrow.hh"
#include "Mu2eG4/inc/constructPS.hh"
#include "Mu2eG4/inc/constructPSShield.hh"
#include "Mu2eG4/inc/nestTubs.hh"
#include "Mu2eG4/inc/finishNesting.hh"

// G4 includes
#include "G4ThreeVector.hh"
#include "G4Material.hh"
#include "G4Color.hh"
#include "G4Polycone.hh"

#include "CLHEP/Units/SystemOfUnits.h"

using namespace std;

namespace mu2e {

  void constructPS(VolumeInfo const & parent, SimpleConfig const & _config) {
    
    // Extract information from the config file.

    GeomHandle<Beamline> beamg;
    double solenoidOffset = beamg->solenoidOffset();
    double rTorus         = beamg->getTS().torusRadius();
    double ts1HalfLength  = beamg->getTS().getTS1().getHalfLength();

    // double ps1HalfLength     = _config.getDouble("toyPS1.vacHalfLength");

    // Build the barrel of the cryostat.
//     TubsParams psCryoParams( _config.getDouble("toyPS.rIn"),
//                              _config.getDouble("toyPS.rOut"),
//                              _config.getDouble("toyPS.CryoHalfLength"));

    // In the Mu2e coordinate system.
    //    double psCryoZ0 = -rTorus + -2.*ts1HalfLength - psCryoParams.zHalfLength();
    // G4ThreeVector psCryoPosition( solenoidOffset, 0., psCryoZ0 );


    double _psLocalOriginZ;
    int    verbosityLevel;
    bool   _psVisible;
    bool   _psSolid;

    double _psVacVesselrIn;
    double _psVacVesselrOut;
    double _psVacVesselWallThickness;
    double _psVacVesselEndPlateHalfThickness;
    double _psVacVesselHalfLength;
    std::string _psVacVesselMaterialName;
    std::string _psInsideMaterialName;

    _psLocalOriginZ                   = _config.getDouble("PS.localOriginZ"); // will not use for now
    // FIXME use the above once params are reconciled

    verbosityLevel                   = _config.getInt("PS.verbosityLevel");
    _psVisible                        = _config.getBool("PS.visible");
    _psSolid                          = _config.getBool("PS.solid");
                                                
    _psVacVesselrIn                   = _config.getDouble("PS.VacVessel.rIn");
    _psVacVesselrOut                  = _config.getDouble("PS.VacVessel.rOut");
    _psVacVesselWallThickness         = _config.getDouble("PS.VacVessel.WallThickness");
    _psVacVesselHalfLength            = _config.getDouble("PS.VacVessel.HalfLength");
    _psVacVesselEndPlateHalfThickness = _config.getDouble("PS.VacVessel.EndPlateHalfThickness");
    _psVacVesselMaterialName          = _config.getString("PS.VacVessel.materialName");
    _psInsideMaterialName             = _config.getString("PS.insideMaterialName");

    G4Material* psVacVesselMaterial = findMaterialOrThrow(_psVacVesselMaterialName);
    G4Material* psInsideMaterial    = findMaterialOrThrow(_psInsideMaterialName);

    double        psCenterZ0 = -rTorus + -2.*ts1HalfLength - _psVacVesselHalfLength;
    G4ThreeVector psMu2eOffset( solenoidOffset, 0., psCenterZ0 );

    double psCenterZ0BasedOnLocalOriginZ = _psLocalOriginZ + 220. - 212 - 30. + _psVacVesselHalfLength;

    verbosityLevel >0 && 
      cout << __func__ << " verbosityLevel                   : " << verbosityLevel  << endl;

    //    bool toyPSVisible        = _config.getBool("toyPS.visible",true);
    //    bool toyPSSolid          = _config.getBool("toyPS.solid",true);
    bool forceAuxEdgeVisible = _config.getBool("g4.forceAuxEdgeVisible",false);
    bool doSurfaceCheck      = _config.getBool("g4.doSurfaceCheck",false);
    //    doSurfaceCheck      = true;
    bool const placePV       = true;

    G4ThreeVector _hallOriginInMu2e = parent.centerInMu2e();

    // we need to replace this in future getting/moving the info from/to geom service

    // VacVessel

    // VacVessel is a "G4Tubs"; or a combination of those

    // Build the barrel of the cryostat, the length is the "outer length"

    TubsParams psVacVesselInnerParams( _psVacVesselrIn,
                                       _psVacVesselrIn +_psVacVesselWallThickness,
                                       _psVacVesselHalfLength);


    VolumeInfo psVacVesselInnerInfo   = nestTubs( "PSVacVesselInner",
                                                  psVacVesselInnerParams,
                                                  psVacVesselMaterial,
                                                  0,
                                                  psMu2eOffset-_hallOriginInMu2e,
                                                  parent,
                                                  0,
                                                  _psVisible,
                                                  G4Colour::Green(),
                                                  _psSolid,
                                                  forceAuxEdgeVisible,
                                                  placePV,
                                                  doSurfaceCheck
                                                  );

    TubsParams psVacVesselOuterParams( _psVacVesselrOut-_psVacVesselWallThickness,
                                       _psVacVesselrOut,
                                       _psVacVesselHalfLength);

    VolumeInfo psVacVesselOuterInfo   = nestTubs( "PSVacVesselOuter",
                                                  psVacVesselOuterParams,
                                                  psVacVesselMaterial,
                                                  0,
                                                  psMu2eOffset-_hallOriginInMu2e,
                                                  parent,
                                                  0,
                                                  _psVisible,
                                                  G4Colour::Green(),
                                                  _psSolid,
                                                  forceAuxEdgeVisible,
                                                  placePV,
                                                  doSurfaceCheck
                                                  );

    // there will be two of them
    TubsParams psVacVesselEndPlateParams( _psVacVesselrIn +_psVacVesselWallThickness,
                                          _psVacVesselrOut-_psVacVesselWallThickness,
                                          _psVacVesselEndPlateHalfThickness);

    double        psVacVesselEndPlateDMu2eOffsetZ0 = 
      psCenterZ0+_psVacVesselHalfLength-_psVacVesselEndPlateHalfThickness;
    G4ThreeVector psVacVesselEndPlateDMu2eOffset( solenoidOffset, 0., psVacVesselEndPlateDMu2eOffsetZ0 );

    VolumeInfo psVacVesselEndPlateDInfo = nestTubs( "PSVacVesselEndPlateD",
                                                    psVacVesselEndPlateParams,
                                                    psVacVesselMaterial,
                                                    0,
                                                    psVacVesselEndPlateDMu2eOffset-_hallOriginInMu2e,
                                                    parent,
                                                    0,
                                                    _psVisible,
                                                    G4Colour::Yellow(),
                                                    _psSolid,
                                                    forceAuxEdgeVisible,
                                                    placePV,
                                                    doSurfaceCheck
                                                    );

    double        psVacVesselEndPlateUMu2eOffsetZ0 = 
      psCenterZ0-_psVacVesselHalfLength+_psVacVesselEndPlateHalfThickness;
    G4ThreeVector psVacVesselEndPlateUMu2eOffset( solenoidOffset, 0., psVacVesselEndPlateUMu2eOffsetZ0 );

    VolumeInfo psVacVesselEndPlateUInfo = nestTubs( "PSVacVesselEndPlateU",
                                                    psVacVesselEndPlateParams,
                                                    psVacVesselMaterial,
                                                    0,
                                                    psVacVesselEndPlateUMu2eOffset-_hallOriginInMu2e,
                                                    parent,
                                                    0,
                                                    _psVisible,
                                                    G4Colour::Red(),
                                                    _psSolid,
                                                    forceAuxEdgeVisible,
                                                    placePV,
                                                    doSurfaceCheck
                                                    );


    //Coil "Outer Shell"

    double _psCoilShell1zOffset;

    double _psCoilShellrIn;
    std::string _psCoilShellMaterialName;

    // Z offset from the local  origin
    double _psCoilShell1zGap;
    // outer radius
    double _psCoilShell1rOut;
    double _psCoilShell1Length;

    // offset from coilShell1
    double _psCoilShell2zGap;
    double _psCoilShell2rOut;
    double _psCoilShell2Length;

    // offset from coilShell2
    double _psCoilShell3zGap;
    double _psCoilShell3rOut;
    double _psCoilShell3Length;

// coil "outer shell
    _psCoilShellrIn                   = _config.getDouble("PS.CoilShell.rIn");
    _psCoilShellMaterialName          = _config.getString("PS.CoilShell.materialName");
                                                
// Z offset from the local  origin
    _psCoilShell1zOffset              = _config.getDouble("PS.CoilShell1.zOffset");
    _psCoilShell1zGap                 = _config.getDouble("PS.CoilShell1.zGap");
// outer radius
    _psCoilShell1rOut                 = _config.getDouble("PS.CoilShell1.rOut");
    _psCoilShell1Length               = _config.getDouble("PS.CoilShell1.Length");
                                                
// offset from coilShell1
    _psCoilShell2zGap                 = _config.getDouble("PS.CoilShell2.zGap");
    _psCoilShell2rOut                 = _config.getDouble("PS.CoilShell2.rOut");
    _psCoilShell2Length               = _config.getDouble("PS.CoilShell2.Length");
                                                
// offset from coilShell2
    _psCoilShell3zGap                 = _config.getDouble("PS.CoilShell3.zGap");
    _psCoilShell3rOut                 = _config.getDouble("PS.CoilShell3.rOut");
    _psCoilShell3Length               = _config.getDouble("PS.CoilShell3.Length");


    //    PolyConeParams psCoilShellParams

    double & ria  = _psCoilShellrIn;

    double  rInner[] = {ria,ria,ria,ria};

    double  rOuter[] = {_psCoilShell1rOut,
                        _psCoilShell1rOut,
                        _psCoilShell2rOut,
                        _psCoilShell3rOut};

    double z1 = _psCoilShell1zGap   + _psCoilShell1Length + _psCoilShell2zGap;
    double z2 = _psCoilShell2Length + _psCoilShell3zGap   + _psCoilShell3Length;

    double  zPlanes[] = {0.,
                         z1,
                         z1,
                         z1+z2};

    unsigned int nPlanes = sizeof(rOuter)/sizeof(rOuter[0]);

    assert ( nPlanes == (sizeof(rOuter) /sizeof( rOuter[0])) );
    assert ( nPlanes == (sizeof(zPlanes)/sizeof(zPlanes[0])) );
    
    double psCoilShellLength = z1+z2;

    //we will place the shell inside the parent, which is the hall
    string const psCoilShellName = "PSCoilShell";
  
    // note that the origin of the G4Polycone is not at its center but
    // at the center of the edge wall of lower Z

    G4ThreeVector psCoilShellLocalOffset = psMu2eOffset - parent.centerInMu2e() -
      G4ThreeVector( 0., 0., psCoilShellLength*0.5 );

    if ( verbosityLevel > 0) {
      cout << __func__ << " solenoidOffset                   : " << solenoidOffset  << endl;
      cout << __func__ << " psCenterZ0                       : " << psCenterZ0      << endl;
      cout << __func__ << " _psLocalOriginZ                  : " << _psLocalOriginZ << endl;
      cout << __func__ << " _psVacVesselHalfLength           : " << _psVacVesselHalfLength << endl;
      cout << __func__ << " psCenterZ0BasedOnLocalOriginZ    : " << psCenterZ0BasedOnLocalOriginZ << endl;
      cout << __func__ << " psCoilShellLength                : " << psCoilShellLength      << endl;
    }

    VolumeInfo psCoilShellInfo(psCoilShellName,psCoilShellLocalOffset,parent.centerInWorld);

    psCoilShellInfo.solid  =  new G4Polycone( psCoilShellName,
                                              0.,
                                              CLHEP::twopi,
                                              nPlanes,
                                              zPlanes,
                                              rInner,
                                              rOuter);
    
    G4Material* psCoilShellMaterial = findMaterialOrThrow(_psCoilShellMaterialName);

    finishNesting(psCoilShellInfo,
                  psCoilShellMaterial,
                  0,
                  psCoilShellLocalOffset,
                  parent.logical,
                  0,
                  _psVisible,
                  G4Colour::White(),
                  _psSolid,
                  forceAuxEdgeVisible,
                  placePV,
                  doSurfaceCheck
                  );

    // the superconducting Coils

    double _psCoilrIn;
    std::string _psCoilMaterialName;

    // Z offset from the local origin
    double _psCoil1zOffset;
    double _psCoil1zGap;
    // outer radius
    double _psCoil1rOut;
    double _psCoil1Length;

    // offset from coil1
    double _psCoil2zGap;
    double _psCoil2rOut;
    double _psCoil2Length;

    // offset from coil2
    double _psCoil3zGap;
    double _psCoil3rOut;
    double _psCoil3Length;

// the three coils are done "together"
// all have the same inner radius
    _psCoilrIn                        = _config.getDouble("PS.Coil.rIn");
    _psCoilMaterialName               = _config.getString("PS.Coil.materialName");
// Z offset from the local origin

    _psCoil1zOffset                   = _config.getDouble("PS.Coil1.zOffset");
    _psCoil1zGap                      = _config.getDouble("PS.Coil1.zGap");
// outer radius
    _psCoil1rOut                      = _config.getDouble("PS.Coil1.rOut");
    _psCoil1Length                    = _config.getDouble("PS.Coil1.Length");
                                                
// offset from coil1
    _psCoil2zGap                      = _config.getDouble("PS.Coil2.zGap");
    _psCoil2rOut                      = _config.getDouble("PS.Coil2.rOut");
    _psCoil2Length                    = _config.getDouble("PS.Coil2.Length");
                                                
// offset from coil2
    _psCoil3zGap                      = _config.getDouble("PS.Coil3.zGap");
    _psCoil3rOut                      = _config.getDouble("PS.Coil3.rOut");
    _psCoil3Length                    = _config.getDouble("PS.Coil3.Length");

    //the coils will be implemented as "G4Tubs" placed inside the shell 

    G4Material* psCoilMaterial = findMaterialOrThrow(_psCoilMaterialName);

    TubsParams psCoil1Params( _psCoilrIn, _psCoil1rOut, _psCoil1Length*0.5);
    
    G4ThreeVector psCoil1LocalOffset(0.,0.,_psCoil1zGap + _psCoil1Length*0.5);

    VolumeInfo psCoil1Info   = nestTubs( "PSCoil1",
                                         psCoil1Params,
                                         psCoilMaterial,
                                         0,
                                         psCoil1LocalOffset,
                                         psCoilShellInfo,
                                         0,
                                         _psVisible,
                                         G4Colour::Red(),
                                         _psSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    TubsParams psCoil2Params( _psCoilrIn, _psCoil2rOut, _psCoil2Length*0.5);
    
    G4ThreeVector psCoil2LocalOffset(0.,0.,
                                     _psCoil1zGap + _psCoil1Length +
                                     _psCoil2zGap + _psCoil2Length*0.5);

    VolumeInfo psCoil2Info   = nestTubs( "PSCoil2",
                                         psCoil2Params,
                                         psCoilMaterial,
                                         0,
                                         psCoil2LocalOffset,
                                         psCoilShellInfo,
                                         0,
                                         _psVisible,
                                         G4Colour::Red(),
                                         _psSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    TubsParams psCoil3Params( _psCoilrIn, _psCoil3rOut, _psCoil3Length*0.5);
    
    G4ThreeVector psCoil3LocalOffset(0.,0.,
                                     _psCoil1zGap + _psCoil1Length +
                                     _psCoil2zGap + _psCoil2Length +
                                     _psCoil3zGap + _psCoil3Length*0.5);

    VolumeInfo psCoil3Info   = nestTubs( "PSCoil3",
                                         psCoil3Params,
                                         psCoilMaterial,
                                         0,
                                         psCoil3LocalOffset,
                                         psCoilShellInfo,
                                         0,
                                         _psVisible,
                                         G4Colour::Red(),
                                         _psSolid,
                                         forceAuxEdgeVisible,
                                         placePV,
                                         doSurfaceCheck
                                         );

    // Build the main PS vacuum body.

    // we may need to adjust it to make sure it conforms to the newer
    // drawings where the collimators enter the PS area

    TubsParams psVacParams( 0.,
                            _psVacVesselrIn,
                            _psVacVesselHalfLength);

    VolumeInfo psVacInfo   = nestTubs( "PSVacuum",
                                       psVacParams,
                                       psInsideMaterial,
                                       0,
                                       psMu2eOffset-_hallOriginInMu2e,
                                       parent,
                                       0,
                                       false,
                                       // _psVisible,
                                       G4Colour::Blue(),
                                       _psSolid,
                                       forceAuxEdgeVisible,
                                       placePV,
                                       doSurfaceCheck
                                       );

    // Build the production target.
    GeomHandle<ProductionTarget> tgt;
    TubsParams prodTargetParams( 0., tgt->rOut(), tgt->halfLength());

    G4Material* prodTargetMaterial = findMaterialOrThrow(_config.getString("targetPS_materialName"));
    
    bool prodTargetVisible = _config.getBool("targetPS.visible",true);
    bool prodTargetSolid   = _config.getBool("targetPS.solid",true);

    VolumeInfo prodTargetInfo   = nestTubs( "ProductionTarget",
                                            prodTargetParams,
                                            prodTargetMaterial,
                                            &tgt->productionTargetRotation(),
                                            tgt->position() - psMu2eOffset,
                                            psVacInfo,
                                            0,
                                            prodTargetVisible,
                                            G4Colour::Magenta(),
                                            prodTargetSolid,
                                            forceAuxEdgeVisible,
                                            placePV,
                                            doSurfaceCheck
                                            );

    constructPSShield(psVacInfo, _config);

    // To compare with g4beamline studies: close the vacuum with a solid disk
    // FIXME
    const double toyEnclosureThickness(_config.getDouble("toyPS.toyEnclosure.Thickness", 0)*CLHEP::mm);
    if(toyEnclosureThickness > 0) {
      
      TubsParams diskParam(0, _config.getDouble("toyPS.rOut"), 0.5*toyEnclosureThickness);
      CLHEP::Hep3Vector toyEnclosurePosition = psMu2eOffset -
        CLHEP::Hep3Vector(0, 0, psVacVesselInnerParams.zHalfLength()) -
        CLHEP::Hep3Vector(0, 0, 0.5*toyEnclosureThickness);

      nestTubs( "PSToyEnclosure",
		diskParam,
		findMaterialOrThrow(_config.getString("toyPS.toyEnclosure.materialName")),
		0,
		toyEnclosurePosition - _hallOriginInMu2e,
		parent,
		0,
		_config.getBool("toyPS.toyEnclosure.visible", true),
		G4Colour::Cyan(),
		_config.getBool("toyPS.toyEnclosure.solid", true),
		forceAuxEdgeVisible,
		placePV,
		doSurfaceCheck
		);
    }

  } // end Mu2eWorld::constructPS
}
