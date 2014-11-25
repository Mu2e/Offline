// Geometry of the detector solenoid
//
// The z-position of the DS is determined based on:
//  - the torus radius of the TS
//  - the half-length of TS5
//  - the half-length of the first DS vacuum volume
//  - the half-length of the front face of the solenoid
//  - the half-length of the DS
//  (see DetectorSolenoidMaker.cc for computation)
//
// Original author Kyle Knoepfel, 2013

#ifndef DETECTORSOLENOID_HH
#define DETECTORSOLENOID_HH

#include <vector>

// #include "art/Persistency/Common/Wrapper.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class DetectorSolenoidMaker;

  class DetectorSolenoid : virtual public Detector {
  public:

    // DS cryostat
    std::string material() const { return _materialName; }
    std::string insideMaterial() const { return _insideMaterialName; }
    double rIn1()  const { return _rIn1; } // inner radius of inner part of DS cryo.
    double rIn2()  const { return _rIn2; } // outer radius of inner part of DS cryo.
    double rOut1() const { return _rOut1; }  // inner radius of outer part of DS cryo.
    double rOut2() const { return _rOut2; }  // outer radius of outer part of DS cryo.
    double halfLength() const { return _halfLength; }
    double endWallHalfLength() const { return _endWallHalfLength; }
    double frontHalfLength() const { return _frontHalfLength; }
    double cryoZMax() const { return _position[CLHEP::Hep3Vector::Z] + _halfLength; }
    const CLHEP::Hep3Vector& position() const { return _position; } // in mu2e coordinates

    // DS shield
    std::string shield_material() const { return _shield_materialName; }
    std::string shield_insideMaterial() const { return _shield_insideMaterialName; }
    double shield_zOffset() const { return _shield_zOffset; } // - this is a superfluous quantity since
    double shield_halfLength() const { return _shield_halfLength; } // the thermal shield is positioned in the same
    double shield_endWallHalfLength() const { return _shield_endWallHalfLength; } // place as the cryo above
    double shield_rIn1() const { return _shield_rIn1; }
    double shield_rIn2() const { return _shield_rIn2; }
    double shield_rOut1() const { return _shield_rOut1; }
    double shield_rOut2() const { return _shield_rOut2; }

    // DS solenoid coils
    std::string coil_material() const { return _coil_materialName; }
    double coil_rIn() const { return _coil_rIn; }
    int nCoils() const { return _nCoils; }
    const std::vector<double>& coil_rOut() const { return _coil_rOut; }
    const std::vector<double>& coil_zLength() const { return _coil_zLength; }
    const std::vector<double>& coil_zPosition() const { return _coil_zPosition; }

    // DS coil support system
    std::string support_material() const { return _support_materialName; }
    double support_rIn() const { return _support_rIn; }
    double support_rOut() const { return _support_rOut; }
    double support_halfLength() const { return _support_halfLength; }

    // Vacuum volumes inside DS
    //
    // The subdivision of the DS vacuum volume is not physical,
    // but it needs to be in Geometry because real physical
    // pieces are placed inside.
    double vac_halfLengthDs1() const { return _ds1HalfLength; }
    double vac_halfLengthDs2() const { return _ds2HalfLength; }
    double vac_zLocDs23Split() const { return _locationDs23Split; }
    double vac_pressure()      const { return _vacuumPressure;    }
    std::string vacuumMaterial() const { return _vacuumMaterialName; }


    //----------------------------------------------------------------
  private:
    friend class DetectorSolenoidMaker;

    // Private ctr: the class should be only obtained via DetectorSolenoid::DetectorSolenoidMaker.
    DetectorSolenoid();

    // DS cryostat
    std::string _materialName; 
    std::string _insideMaterialName; 
    double _rIn1; 
    double _rIn2; 
    double _rOut1;
    double _rOut2;
    double _halfLength; 
    double _endWallHalfLength; 
    double _frontHalfLength; 
    CLHEP::Hep3Vector _position; 

    // DS shield
    std::string _shield_materialName;
    std::string _shield_insideMaterialName;
    double _shield_zOffset;
    double _shield_halfLength;
    double _shield_endWallHalfLength;
    double _shield_rIn1;
    double _shield_rIn2;
    double _shield_rOut1;
    double _shield_rOut2;

    // DS solenoid coils
    const int _nCoils = 11;
    std::string         _coil_materialName; 
    double              _coil_rIn; 
    std::vector<double> _coil_rOut; 
    std::vector<double> _coil_zLength; 
    std::vector<double> _coil_zPosition; 

    // DS coil support system
    std::string _support_materialName; 
    double      _support_rIn; 
    double      _support_rOut; 
    double      _support_halfLength; 

    // Vacuum volumes inside DS
    double _vacuumPressure;
    std::string _vacuumMaterialName; 
    double _ds1HalfLength;
    double _ds2HalfLength;
    double _locationDs23Split;

    // Needed for persistency
    //    template<class T> friend class art::Wrapper;
    //    DetectorSolenoid() {}
  };
}

#endif/*DETECTORSOLENOID_HH*/
