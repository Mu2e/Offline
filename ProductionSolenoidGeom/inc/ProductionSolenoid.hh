#ifndef ProductionSolenoidGeom_ProductionSolenoid_hh
#define ProductionSolenoidGeom_ProductionSolenoid_hh


//
//
// Original author KLG
//
// Added code for rings at saddle points

#include "CLHEP/Vector/ThreeVector.h"
#include "Offline/GeomPrimitives/inc/Tube.hh"
#include "Offline/GeomPrimitives/inc/Polycone.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class ProductionSolenoidMaker;

  class ProductionSolenoid : virtual public Detector {

  public:
    ~ProductionSolenoid() override = default;

    // delete  automatic copy/assignments as not needed (would be incorrect due to unique_ptr anyway)
    ProductionSolenoid( ProductionSolenoid const & ) = delete;
    ProductionSolenoid( ProductionSolenoid&&       ) = delete;
    ProductionSolenoid& operator=( ProductionSolenoid const & ) = delete;
    ProductionSolenoid& operator=( ProductionSolenoid&&       ) = delete;

    // do we need more than that? does the Tube have all the accessors?

    Tube const * getVacVesselInnerParamsPtr() const {return _psVacVesselInnerParams.get();}
    Tube const * getVacVesselOuterParamsPtr() const {return _psVacVesselOuterParams.get();}

    Tube const * getVacVesselEndPlateDParamsPtr() const { return _psVacVesselEndPlateDParams.get(); }
    Tube const * getVacVesselEndPlateUParamsPtr() const { return _psVacVesselEndPlateUParams.get(); }

    Polycone const * getCoilShellParamsPtr() const { return _psCoilShellParams.get(); }

    Tube const * getCoil1ParamsPtr() const { return _psCoil1Params.get(); }
    Tube const * getCoil2ParamsPtr() const { return _psCoil2Params.get(); }
    Tube const * getCoil3ParamsPtr() const { return _psCoil3Params.get(); }

    Tube const * getRing1ParamsPtr() const { return _psRing1Params.get(); }
    Tube const * getRing2ParamsPtr() const { return _psRing2Params.get(); }

    // The point on the PS axis at the end of the cryostat on the downstream (proton exit) side
    const CLHEP::Hep3Vector& psEndRefPoint() const { return _psEndRefPoint; }

  private:

    friend class ProductionSolenoidMaker;

    // The class should only be constructed via ProductionSolenoid::ProductionSolenoidMaker.
    ProductionSolenoid(){};

    // it has several components

    // VacVessel is a set of Tubes; lets try to put all of it here

    std::unique_ptr<Tube> _psVacVesselInnerParams;
    std::unique_ptr<Tube> _psVacVesselOuterParams;

    std::unique_ptr<Tube> _psVacVesselEndPlateDParams;
    std::unique_ptr<Tube> _psVacVesselEndPlateUParams;

    // The rings are tubes
    std::unique_ptr<Tube> _psRing1Params;
    std::unique_ptr<Tube> _psRing2Params;

    // CoilShell is a Polycone
    std::unique_ptr<Polycone> _psCoilShellParams;

    // Coils are Tubes

    std::unique_ptr<Tube> _psCoil1Params;
    std::unique_ptr<Tube> _psCoil2Params;
    std::unique_ptr<Tube> _psCoil3Params;

    CLHEP::Hep3Vector   _psEndRefPoint;
  };
}

#endif/*ProductionSolenoidGeom_ProductionSolenoid_hh*/
