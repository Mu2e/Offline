#ifndef ProductionSolenoidGeom_ProductionSolenoid_hh
#define ProductionSolenoidGeom_ProductionSolenoid_hh


//
// $Id: ProductionSolenoid.hh,v 1.1 2012/03/13 19:04:11 genser Exp $
// $Author: genser $
// $Date: 2012/03/13 19:04:11 $
//
// Original author KLG
//

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "GeomPrimitives/inc/Tube.hh"
#include "GeomPrimitives/inc/Polycone.hh"
#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {
    
  class ProductionSolenoidMaker;

  class ProductionSolenoid : virtual public Detector {

  public: 

    // do we need more than that? does the Tube have all the accessors?

    Tube const * getVacVesselInnerParamsPtr() const {return _psVacVesselInnerParams.get();}
    Tube const * getVacVesselOuterParamsPtr() const {return _psVacVesselOuterParams.get();}

    Tube const * getVacVesselEndPlateDParamsPtr() const { return _psVacVesselEndPlateDParams.get(); }
    Tube const * getVacVesselEndPlateUParamsPtr() const { return _psVacVesselEndPlateUParams.get(); }

    Polycone const * getCoilShellParamsPtr() const { return _psCoilShellParams.get(); }

    Tube const * getCoil1ParamsPtr() const { return _psCoil1Params.get(); }
    Tube const * getCoil2ParamsPtr() const { return _psCoil2Params.get(); }
    Tube const * getCoil3ParamsPtr() const { return _psCoil3Params.get(); }

    Tube const * getVacuumParamsPtr() const { return _psVacuumParams.get(); }

    Tube const * getToyEnclosureParamsPtr() const { return _psToyEnclosureParams.get(); }

  private: 

    friend class ProductionSolenoidMaker;

    // The class should only be constructed via ProductionSolenoid::ProductionSolenoidMaker.
    ProductionSolenoid(){};

    // hide automatic copy/assignments as not needed (would be incorrect due to auto_ptr anyway)
    ProductionSolenoid( ProductionSolenoid const & );
    ProductionSolenoid const & operator= ( ProductionSolenoid const & );

    // it has several components

    // VacVessel is a set of Tubes; lets try to put all of it here

    std::auto_ptr<Tube> _psVacVesselInnerParams;
    std::auto_ptr<Tube> _psVacVesselOuterParams;

    std::auto_ptr<Tube> _psVacVesselEndPlateDParams;
    std::auto_ptr<Tube> _psVacVesselEndPlateUParams;

    // CoilShell is a Polycone
    std::auto_ptr<Polycone> _psCoilShellParams;

    // Coils are Tubes

    std::auto_ptr<Tube> _psCoil1Params;
    std::auto_ptr<Tube> _psCoil2Params;
    std::auto_ptr<Tube> _psCoil3Params;

    // Vacuum
    std::auto_ptr<Tube> _psVacuumParams;

    // ToyEnclosure
    std::auto_ptr<Tube> _psToyEnclosureParams;

  };
}

#endif/*ProductionSolenoidGeom_ProductionSolenoid_hh*/
