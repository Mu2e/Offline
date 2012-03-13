#ifndef ProductionSolenoidGeom_ProductionSolenoidMaker_hh
#define ProductionSolenoidGeom_ProductionSolenoidMaker_hh
//
// Class to construct and return ProductionSolenoid
//
// $Id: ProductionSolenoidMaker.hh,v 1.1 2012/03/13 19:04:11 genser Exp $
// $Author: genser $
// $Date: 2012/03/13 19:04:11 $
//
// Original author KLG
//

#include <memory>
#include <string>

#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class ProductionSolenoid;
  class SimpleConfig;
  class Tube;

  class ProductionSolenoidMaker {

  public:

    ProductionSolenoidMaker( SimpleConfig const & config, 
                             double solenoidOffset,
                             double rTorus,
                             double ts1HalfLength);

    void parseConfig( SimpleConfig const & _config );

    // This is deprecated and will go away soon.
    // Still needed for root graphics version.
    const ProductionSolenoid& getProductionSolenoid() const { return *_ps;}

    // This is the accessor that will remain.
    std::auto_ptr<ProductionSolenoid> getProductionSolenoidPtr() { return _ps; }

  private:

    // hide automatic copy/assignments as not needed (would be incorrect due to auto_ptr anyway)
    ProductionSolenoidMaker( ProductionSolenoidMaker const & );
    ProductionSolenoidMaker const & operator= ( ProductionSolenoidMaker const & );

    std::auto_ptr<ProductionSolenoid> _ps;

    double _psLocalOriginZ;
    int    _verbosityLevel;
    bool   _psVisible;
    bool   _psSolid;

    double _psVacVesselrIn;
    double _psVacVesselrOut;
    double _psVacVesselWallThickness;
    double _psVacVesselEndPlateHalfThickness;
    double _psVacVesselHalfLength;
    std::string _psVacVesselMaterialName;
    std::string _psInsideMaterialName;

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

    double _psToyEnclosureThickness;
    std::string _psToyEnclosureMaterialName;

    // derived quantities etc...

  };

}  //namespace mu2e

#endif /* ProductionSolenoidGeom_ProductionSolenoidMaker_hh */
