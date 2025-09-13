#ifndef STMGeom_SSCSupport_hh
#define STMGeom_SSCSupport_hh

// STM Collimator Object
//
// Author: Haichuan Cao
// Sept 2023

#include <string>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class SSCSupport {
  public:

    SSCSupport(bool build,
                 double table_L, double table_H, double table_T,
                 double leg_L, double leg_H, double leg_T,
                 double base_L, double base_H, double base_T,
                 double wall_L, double wall_H, double wall_T,
                 double hole_H, double hole_T,
                 double FLeadStand_L, double FLeadStand_H, double FLeadStand_T,
                 double FLeadShim_H, double FLeadShim_T,
                 double FAluminumShim_T, double FAluminumExtra_L, double FAluminumExtra_H,
                 CLHEP::Hep3Vector const & originInMu2e = CLHEP::Hep3Vector(),
                 CLHEP::HepRotation const & rotation = CLHEP::HepRotation()
                 ) :
      _build(build),
      _table_L(table_L), _table_H(table_H), _table_T(table_T),
      _leg_L(leg_L), _leg_H(leg_H), _leg_T(leg_T),
      _base_L(base_L), _base_H(base_H), _base_T(base_T),
      _wall_L(wall_L), _wall_H(wall_H), _wall_T(wall_T),
      _hole_H(hole_H), _hole_T(hole_T),
      _FLeadStand_L(FLeadStand_L), _FLeadStand_H(FLeadStand_H), _FLeadStand_T(FLeadStand_T),
      _FLeadShim_H(FLeadShim_H), _FLeadShim_T(FLeadShim_T),
      _FAluminumShim_T(FAluminumShim_T), _FAluminumExtra_L(FAluminumExtra_L), _FAluminumExtra_H(FAluminumExtra_H),
      _originInMu2e(originInMu2e),
      _rotation(rotation)
    {
    }

    bool   build()       const {return _build;}

    double table_L()   const {return _table_L;}
    double table_H()   const {return _table_H;}
    double table_T()   const {return _table_T;}
    double leg_L()   const {return _leg_L;}
    double leg_H()   const {return _leg_H;}
    double leg_T()   const {return _leg_T;}
    double base_L()   const {return _base_L;}
    double base_H()   const {return _base_H;}
    double base_T()   const {return _base_T;}
    double wall_L()   const {return _wall_L;}
    double wall_H()   const {return _wall_H;}
    double wall_T()   const {return _wall_T;}
    double hole_H()   const {return _hole_H;}
    double hole_T()   const {return _hole_T;}

    double FLeadStand_L() const {return _FLeadStand_L;}
    double FLeadStand_H() const {return _FLeadStand_H;}
    double FLeadStand_T() const {return _FLeadStand_T;}
    double FLeadShim_H() const {return _FLeadShim_H;}
    double FLeadShim_T() const {return _FLeadShim_T;}
    double FAluminumShim_T() const {return _FAluminumShim_T;}
    double FAluminumExtra_L() const {return _FAluminumExtra_L;}
    double FAluminumExtra_H() const {return _FAluminumExtra_H;}

    CLHEP::Hep3Vector const &  originInMu2e()     const { return _originInMu2e; }
    CLHEP::HepRotation const & rotation()         const { return _rotation; }
    // Genreflex can't do persistency of vector<SSCSupport> without a default constructor
    SSCSupport() {}

  private:

    bool   _build;
    bool   _VDbuild;

    double _table_L;
    double _table_H;
    double _table_T;
    double _leg_L;
    double _leg_H;
    double _leg_T;
    double _base_L;
    double _base_H;
    double _base_T;
    double _wall_L;
    double _wall_H;
    double _wall_T;
    double _hole_H;
    double _hole_T;

    double _FLeadStand_L;
    double _FLeadStand_H;
    double _FLeadStand_T;
    double _FLeadShim_H;
    double _FLeadShim_T;
    double _FAluminumShim_T;
    double _FAluminumExtra_L;
    double _FAluminumExtra_H;


    CLHEP::Hep3Vector  _originInMu2e;
    CLHEP::HepRotation _rotation; // wrt to parent volume
  };

}

#endif/*STMGeom_SSCSupport_hh*/
