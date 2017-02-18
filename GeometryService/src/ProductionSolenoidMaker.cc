//
// Construct and return ProductionSolenoid
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
#include "GeometryService/inc/ProductionSolenoidMaker.hh"
#include "ProductionSolenoidGeom/inc/ProductionSolenoid.hh"

#include "ConfigTools/inc/SimpleConfig.hh"

using namespace std;

namespace mu2e {

  // Constructor that gets information from the config file instead of
  // from arguments.
  ProductionSolenoidMaker::ProductionSolenoidMaker(SimpleConfig const & _config,
                                                   double solenoidOffset)
  {

    // if( ! _config.getBool("hasProductionSolenoid",false) ) return;

    // create an empty PS
    _ps = unique_ptr<ProductionSolenoid>(new ProductionSolenoid());

    // access its object through a reference

    ProductionSolenoid & ps = *_ps.get();

    parseConfig(_config);

    ps._psEndRefPoint = CLHEP::Hep3Vector(solenoidOffset, 0, _config.getDouble("PS.cryostatRefZ"));

    CLHEP::Hep3Vector psCryoMu2eOffset = ps._psEndRefPoint + CLHEP::Hep3Vector(0, 0, _psVacVesselHalfLength );
    CLHEP::Hep3Vector psCryoInnerExtOffset(0, 0, -_psVacVesselInHalfExt);

    const double psCoilRefZ = _config.getDouble("PS.coilRefZ");

    // now create the specific components

    // Build the barrel of the cryostat, the length is the "outer length"

    ps._psVacVesselInnerParams = std::unique_ptr<Tube>
      (new Tube(_psVacVesselMaterialName,
                psCryoMu2eOffset+psCryoInnerExtOffset,
                 _psVacVesselrIn,
                _psVacVesselrIn +_psVacVesselWallThickness,
                _psVacVesselHalfLength+_psVacVesselInHalfExt ));

    ps._psVacVesselOuterParams = std::unique_ptr<Tube>
      (new Tube(_psVacVesselMaterialName,
                psCryoMu2eOffset,
                _psVacVesselrOut-_psVacVesselWallThickness,
                _psVacVesselrOut,
                _psVacVesselHalfLength));

    CLHEP::Hep3Vector _psRing1Mu2eOffset(solenoidOffset,0,_psRing1Mu2eOffsetZ);
    ps._psRing1Params = std::unique_ptr<Tube>
      (new Tube(_psRingMaterialName,
                _psRing1Mu2eOffset,
                _psRingIR,
                _psRingOR,
                _psRingLength/2.0));

    CLHEP::Hep3Vector _psRing2Mu2eOffset(solenoidOffset,0,_psRing2Mu2eOffsetZ);
    ps._psRing2Params = std::unique_ptr<Tube>
      (new Tube(_psRingMaterialName,
                _psRing2Mu2eOffset,
                _psRingIR,
                _psRingOR,
                _psRingLength/2.0));

    // two endplates

    CLHEP::Hep3Vector endPlateShift(0, 0, _psVacVesselHalfLength-_psVacVesselEndPlateHalfThickness);

    ps._psVacVesselEndPlateDParams = std::unique_ptr<Tube>
      (new Tube(_psVacVesselMaterialName,
                psCryoMu2eOffset + endPlateShift,
                _psVacVesselrIn +_psVacVesselWallThickness,
                _psVacVesselrOut-_psVacVesselWallThickness,
                _psVacVesselEndPlateHalfThickness));

    ps._psVacVesselEndPlateUParams = std::unique_ptr<Tube>
      (new Tube(_psVacVesselMaterialName,
                psCryoMu2eOffset - endPlateShift,
                _psVacVesselrIn +_psVacVesselWallThickness,
                _psVacVesselrOut-_psVacVesselWallThickness,
                _psVacVesselEndPlateHalfThickness));


    //    PolyConeParams psCoilShellParams
    std::vector<double> zPlane, rIn, rOut;
    _config.getVectorDouble("PS.CoilShell.zPlane", zPlane);
    _config.getVectorDouble("PS.CoilShell.rIn",    rIn, zPlane.size());
    _config.getVectorDouble("PS.CoilShell.rOut",   rOut, zPlane.size());

    //aka originInMu2e:
    CLHEP::Hep3Vector psCoilShellMu2eOffset(solenoidOffset, 0, psCoilRefZ);

    ps._psCoilShellParams = std::unique_ptr<Polycone>
      (new Polycone(zPlane,
                    rIn,
                    rOut,
                    psCoilShellMu2eOffset,
                    _config.getString("PS.CoilShell.materialName")
                    ));

    // Coils are Tubes placed inside the Coil Shell
    CLHEP::Hep3Vector psCoil1Mu2eOffset(solenoidOffset,
                                        0,
                                        psCoilRefZ + _psCoil1zOffset +
                                        _psCoil1Length*0.5);

    ps._psCoil1Params = std::unique_ptr<Tube>
      (new Tube(_psCoilMaterialName,
                psCoil1Mu2eOffset,
                _psCoilrIn, _psCoil1rOut, _psCoil1Length*0.5));

    CLHEP::Hep3Vector psCoil2Mu2eOffset(solenoidOffset,
                                        0,
                                        psCoilRefZ + _psCoil1zOffset +
                                        _psCoil1Length + _psCoil2zGap +
                                        _psCoil2Length*0.5);

    ps._psCoil2Params = std::unique_ptr<Tube>
      (new Tube(_psCoilMaterialName,
                psCoil2Mu2eOffset,
                _psCoilrIn, _psCoil2rOut, _psCoil2Length*0.5));

    CLHEP::Hep3Vector psCoil3Mu2eOffset(solenoidOffset,
                                        0,
                                        psCoilRefZ + _psCoil1zOffset +
                                        _psCoil1Length + _psCoil2zGap +
                                        _psCoil2Length + _psCoil3zGap +
                                        _psCoil3Length*0.5);

    ps._psCoil3Params = std::unique_ptr<Tube>
      (new Tube(_psCoilMaterialName,
                psCoil3Mu2eOffset,
                _psCoilrIn, _psCoil3rOut, _psCoil3Length*0.5));

  }

  void ProductionSolenoidMaker::parseConfig( SimpleConfig const & _config ){

    _verbosityLevel                   = _config.getInt("PS.verbosityLevel");

    _psVacVesselrIn                   = _config.getDouble("PS.VacVessel.rIn");
    _psVacVesselrOut                  = _config.getDouble("PS.VacVessel.rOut");
    _psVacVesselWallThickness         = _config.getDouble("PS.VacVessel.WallThickness");
    _psVacVesselHalfLength            = _config.getDouble("PS.VacVessel.HalfLength");
    _psVacVesselEndPlateHalfThickness = _config.getDouble("PS.VacVessel.EndPlateHalfThickness");
    _psVacVesselMaterialName          = _config.getString("PS.VacVessel.materialName");
    _psVacVesselInHalfExt             = _config.getDouble("PS.VacVessel.InnerExtension",0.0)/2.0;

    // Rings added by David Norvil Brown, March 2015
    _psRingMaterialName = _config.getString("PS.ring.materialName");
    _psRingOR = _config.getDouble("PS.ring.rOut");
    _psRingIR = _config.getDouble("PS.ring.rIn");
    _psRingLength = _config.getDouble("PS.ring.length");
    _psRing1Mu2eOffsetZ = _config.getDouble("PS.ring.Z1");
    _psRing2Mu2eOffsetZ = _config.getDouble("PS.ring.Z2");

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
