//
// Construct and return ProductionSolenoid
//
// $Id: ProductionSolenoidMaker.cc,v 1.4 2012/03/29 19:06:31 gandr Exp $
// $Author: gandr $
// $Date: 2012/03/29 19:06:31 $
//
// Original author KLG
//
// Notes

// c++ includes
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

// clhep includes
#include "CLHEP/Vector/ThreeVector.h"

// Framework includes
#include "messagefacility/MessageLogger/MessageLogger.h"

// Mu2e includes
#include "ProductionSolenoidGeom/inc/ProductionSolenoidMaker.hh"
#include "ProductionSolenoidGeom/inc/ProductionSolenoid.hh"

#include "Mu2eUtilities/inc/SimpleConfig.hh"

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  ProductionSolenoidMaker::ProductionSolenoidMaker(SimpleConfig const & _config,
                                                   double solenoidOffset,
                                                   double rTorus,
                                                   double ts1HalfLength)
  {

    // if( ! _config.getBool("hasProductionSolenoid",false) ) return;

    // create an empty PS
    _ps = auto_ptr<ProductionSolenoid>(new ProductionSolenoid());

    // access its object through a reference

    ProductionSolenoid & ps = *_ps.get();

    parseConfig(_config);

    ps._psEndRefPoint = CLHEP::Hep3Vector(solenoidOffset, 0, _config.getDouble("PS.cryostatRefZ"));

    CLHEP::Hep3Vector psCryoMu2eOffset = ps._psEndRefPoint + CLHEP::Hep3Vector(0, 0, _psVacVesselHalfLength );

    const double psCoilRefZ = _config.getDouble("PS.coilRefZ");

    // now create the specific components

    // Build the barrel of the cryostat, the length is the "outer length"

    ps._psVacVesselInnerParams = std::auto_ptr<Tube>
      (new Tube(_psVacVesselMaterialName,
                psCryoMu2eOffset,
                 _psVacVesselrIn,
                _psVacVesselrIn +_psVacVesselWallThickness,
                _psVacVesselHalfLength));

    ps._psVacVesselOuterParams = std::auto_ptr<Tube>
      (new Tube(_psVacVesselMaterialName,
                psCryoMu2eOffset,
                _psVacVesselrOut-_psVacVesselWallThickness,
                _psVacVesselrOut,
                _psVacVesselHalfLength));


    // two endplates

    CLHEP::Hep3Vector endPlateShift(0, 0, _psVacVesselHalfLength-_psVacVesselEndPlateHalfThickness);

    ps._psVacVesselEndPlateDParams = std::auto_ptr<Tube>
      (new Tube(_psVacVesselMaterialName,
                psCryoMu2eOffset + endPlateShift,
                _psVacVesselrIn +_psVacVesselWallThickness,
                _psVacVesselrOut-_psVacVesselWallThickness,
                _psVacVesselEndPlateHalfThickness));

    ps._psVacVesselEndPlateUParams = std::auto_ptr<Tube>
      (new Tube(_psVacVesselMaterialName,
                psCryoMu2eOffset - endPlateShift,
                _psVacVesselrIn +_psVacVesselWallThickness,
                _psVacVesselrOut-_psVacVesselWallThickness,
                _psVacVesselEndPlateHalfThickness));


    //    PolyConeParams psCoilShellParams

    double & ria  = _psCoilShellrIn;

    double  rInner[] = {ria,ria,ria,ria};

    double  rOuter[] = {_psCoilShell1rOut,
                        _psCoilShell1rOut,
                        _psCoilShell2rOut,
                        _psCoilShell3rOut};

    // FIXME: just use polycone params from config
    double z1 = _psCoilShell1zGap   + _psCoilShell1Length + _psCoilShell2zGap;
    double z2 = _psCoilShell2Length + _psCoilShell3zGap   + _psCoilShell3Length;

    double  zPlanes[] = {0.,
                         z1,
                         z1,
                         z1+z2};

    unsigned int nPlanes = sizeof(rOuter)/sizeof(rOuter[0]);

    assert ( nPlanes == (sizeof(rOuter) /sizeof( rOuter[0])) );
    assert ( nPlanes == (sizeof(zPlanes)/sizeof(zPlanes[0])) );

    // Use the cold mass reference point to position the "shell" (stabilizer)
    // Note that the origin of the G4Polycone is not at its center but
    // at the center of the edge wall of lower Z

    //aka originInMu2e:
    CLHEP::Hep3Vector psCoilShellMu2eOffset(solenoidOffset, 0, psCoilRefZ + _psCoil1zOffset - _psCoilShell1zGap );

    // CoilShell is a Polycone

    ps._psCoilShellParams = std::auto_ptr<Polycone>
      (new Polycone(0.,
                    CLHEP::twopi,
                    nPlanes,
                    zPlanes,
                    rInner,
                    rOuter,
                    psCoilShellMu2eOffset,
                    _psCoilShellMaterialName));

    // Coils are Tubes placed inside the Coil Shell
    CLHEP::Hep3Vector psCoil1Mu2eOffset(solenoidOffset,
                                        0,
                                        psCoilRefZ + _psCoil1zOffset +
                                        _psCoil1Length*0.5);

    ps._psCoil1Params = std::auto_ptr<Tube>
      (new Tube(_psCoilMaterialName,
                psCoil1Mu2eOffset,
                _psCoilrIn, _psCoil1rOut, _psCoil1Length*0.5));

    CLHEP::Hep3Vector psCoil2Mu2eOffset(solenoidOffset,
                                        0,
                                        psCoilRefZ + _psCoil1zOffset +
                                        _psCoil1Length + _psCoil2zGap +
                                        _psCoil2Length*0.5);

    ps._psCoil2Params = std::auto_ptr<Tube>
      (new Tube(_psCoilMaterialName,
                psCoil2Mu2eOffset,
                _psCoilrIn, _psCoil2rOut, _psCoil2Length*0.5));

    CLHEP::Hep3Vector psCoil3Mu2eOffset(solenoidOffset,
                                        0,
                                        psCoilRefZ + _psCoil1zOffset +
                                        _psCoil1Length + _psCoil2zGap +
                                        _psCoil2Length + _psCoil3zGap +
                                        _psCoil3Length*0.5);

    ps._psCoil3Params = std::auto_ptr<Tube>
      (new Tube(_psCoilMaterialName,
                psCoil3Mu2eOffset,
                _psCoilrIn, _psCoil3rOut, _psCoil3Length*0.5));

    // PSVacuum
    // we shorten/shift the vacuum to place the vd.  ??? FIXME: why not put vd in the vacuum instead?
    const double vacStartZ =
      ps._psEndRefPoint[2] + (_config.getBool("hasVirtualDetector",false) ? 2*_config.getDouble("vd.halfLength") : 0);

    const double vacEndZ = -rTorus + -2.*ts1HalfLength;

    CLHEP::Hep3Vector psVacuumMu2eOffset(solenoidOffset, 0, 0.5*(vacStartZ + vacEndZ));

    ps._psVacuumParams = std::auto_ptr<Tube>
      (new Tube(_psInsideMaterialName,
                psVacuumMu2eOffset,
                0., _psVacVesselrIn, 0.5*(vacEndZ - vacStartZ)));

  }

  void ProductionSolenoidMaker::parseConfig( SimpleConfig const & _config ){

    _verbosityLevel                   = _config.getInt("PS.verbosityLevel");
    _psVisible                        = _config.getBool("PS.visible");
    _psSolid                          = _config.getBool("PS.solid");

    _psVacVesselrIn                   = _config.getDouble("PS.VacVessel.rIn");
    _psVacVesselrOut                  = _config.getDouble("PS.VacVessel.rOut");
    _psVacVesselWallThickness         = _config.getDouble("PS.VacVessel.WallThickness");
    _psVacVesselHalfLength            = _config.getDouble("PS.VacVessel.HalfLength");
    _psVacVesselEndPlateHalfThickness = _config.getDouble("PS.VacVessel.EndPlateHalfThickness");
    _psVacVesselMaterialName          = _config.getString("PS.VacVessel.materialName");
    _psInsideMaterialName             = _config.getString("PS.insideMaterialName");

    //

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

    //

// the three coils are done "together"
// all have the same inner radius
    _psCoilrIn                        = _config.getDouble("PS.Coil.rIn");
    _psCoilMaterialName               = _config.getString("PS.Coil.materialName");
// Z offset from the local origin

    _psCoil1zOffset                   = _config.getDouble("PS.Coil1.zOffset");
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
  }

} // namespace mu2e
