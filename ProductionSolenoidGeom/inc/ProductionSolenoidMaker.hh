#ifndef ProductionSolenoidGeom_ProductionSolenoidMaker_hh
#define ProductionSolenoidGeom_ProductionSolenoidMaker_hh
//
// Class to construct and return ProductionSolenoid
//
// $Id: ProductionSolenoidMaker.hh,v 1.6 2012/06/06 19:29:31 gandr Exp $
// $Author: gandr $
// $Date: 2012/06/06 19:29:31 $
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
                             double solenoidOffset);

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

    int    _verbosityLevel;

    double _psVacVesselrIn;
    double _psVacVesselrOut;
    double _psVacVesselWallThickness;
    double _psVacVesselEndPlateHalfThickness;
    double _psVacVesselHalfLength;
    std::string _psVacVesselMaterialName;

    // the superconducting Coils

    double _psCoilrIn;
    std::string _psCoilMaterialName;

    // Z offset from the local origin
    double _psCoil1zOffset;
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
