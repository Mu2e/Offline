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
    int coilVersion() const { return _coilVersion; }
    // This is for coilVersion 2
    const std::vector<std::string>& coil_materials() const { return _coil_mats;}
    const std::vector<double>& coil_rOut() const { return _coil_rOut; }
    const std::vector<double>& coil_zLength() const { return _coil_zLength; }
    const std::vector<double>& coil_zPosition() const { return _coil_zPosition; }

    // DS coil spacers
    int nSpacers() const { return _nSpacers; }
    double spacer_rIn() const { return _spacer_rIn; }
    std::string spacer_material() const { return _spacer_materialName; }
    const std::vector<double>& spacer_rOut() const { return _spacer_rOut; }
    const std::vector<double>& spacer_zLength() const { return _spacer_zLength;}
    const std::vector<double>& spacer_zPosition() const { return _spacer_zPosition; }


    // DS coil support system
    std::string support_material() const { return _support_materialName; }
    double support_rIn() const { return _support_rIn; }
    double support_rOut() const { return _support_rOut; }
    double support_halfLength() const { return _support_halfLength; }

    // Support rings - matched to saddles (David Norvil Brown, May 2015)
    double rInRingSide() const { return _rInRingSide; }
    double rOutRingSide() const { return _rOutRingSide; }
    double thickRingSide() const { return _thickRingSide; }
    double rInRing() const { return _rInRing; }
    double rOutRing() const { return _rOutRing; }
    double lengthRing() const { return _lengthRing; }
    std::string RingMaterial() const { return _RingMaterial; }
    std::vector<double> xRing() const { return _xRing; }
    std::vector<double> yRing() const { return _yRing; }
    std::vector<double> zRing() const { return _zRing; }

    // Rails for DS elements to ride on within cryostat
    std::vector<double>  uOutlineRail() const { return _uOutlineRail; }
    std::vector<double>  vOutlineRail() const { return _vOutlineRail; }
    std::string RailMaterial() const { return _RailMaterial; }
    // 2 below refers to part of rail in DS2Vacuum, 3 to DS3Vacuum
    double lengthRail2() const { return _lengthRail2; }
    double lengthRail3() const { return _lengthRail3; }
    CLHEP::Hep3Vector n2RailCenter() const { return _n2RailCenter; }
    CLHEP::Hep3Vector s2RailCenter() const { return _s2RailCenter; }
    CLHEP::Hep3Vector n3RailCenter() const { return _n3RailCenter; }
    CLHEP::Hep3Vector s3RailCenter() const { return _s3RailCenter; }
    // Bearing Blocks for rails, similarly split into DS2Vacuum part and 
    // DS3Vacuum part, except outlines and material  D. No. Brown, Jan 2016.
    std::vector<double>  uOutlineBBlock() const { return _uOutlineBBlock; }
    std::vector<double>  vOutlineBBlock() const { return _vOutlineBBlock; }
    std::string BBlockMaterial() const { return _BBlockMaterial; }
    double lengthBBlock2() const { return _lengthBBlock2; }
    double lengthBBlock3() const { return _lengthBBlock3; }
    std::vector<CLHEP::Hep3Vector> BBlockCenters2() const {
      return _BBlockCenters2; }
    std::vector<CLHEP::Hep3Vector> BBlockCenters3() const {
      return _BBlockCenters3; }

    // MBS Spherical Support outline
    bool hasMBSS() const { return _hasMBSS; }
    double  lengthMBSS()               const { return _lengthMBSS; }
    std::vector<double> uOutlineMBSS() const { return _uOutlineMBSS; }
    std::vector<double> vOutlineMBSS() const { return _vOutlineMBSS; }
    CLHEP::Hep3Vector MBSSlocation()   const { return _locationMBSS; }
    std::string MBSSmaterial()         const { return _materialMBSS; }

    // Cable Runs
    bool hasCableRunCal()              const { return _hasCableRunCal; }
    std::string cableMaterial()        const { return _materialCableRunCal; }
    double rInCableRunCal()            const { return _rInCableRunCal; }
    double rOutCableRunCal()           const { return _rOutCableRunCal; }
    double lengthCableRunCal()         const { return _lengthCableRunCal; }
    double phi0CableRunCal()           const { return _phi0CableRunCal; }
    double dPhiCableRunCal()           const { return _dPhiCableRunCal; }
    double zCCableRunCal()             const { return _zCCableRunCal; }

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
    std::string              _coil_materialName; 
    double                   _coil_rIn; 
    int                      _coilVersion;
    // Allow coil materials to vary in coilVersion 2
    std::vector<std::string> _coil_mats;
    std::vector<double>      _coil_rOut; 
    std::vector<double>      _coil_zLength; 
    std::vector<double>      _coil_zPosition;
 

    // DS coil spacers
    const int                _nSpacers = 5;
    std::string              _spacer_materialName;
    double                   _spacer_rIn;
    std::vector<double>      _spacer_rOut;
    std::vector<double>      _spacer_zLength;
    std::vector<double>      _spacer_zPosition;


    // DS coil support system
    std::string _support_materialName; 
    double      _support_rIn; 
    double      _support_rOut; 
    double      _support_halfLength; 

    // Rings 
    double _rInRingSide, _rOutRingSide, _thickRingSide;
    double _rInRing, _rOutRing, _lengthRing;
    std::string         _RingMaterial;
    std::vector<double> _xRing;
    std::vector<double> _yRing;
    std::vector<double> _zRing;

    // Rails
    std::vector<double>                _uOutlineRail;
    std::vector<double>                _vOutlineRail;
    std::string                        _RailMaterial;
    double                             _lengthRail2;
    double                             _lengthRail3;
    CLHEP::Hep3Vector                  _n2RailCenter;
    CLHEP::Hep3Vector                  _s2RailCenter;
    CLHEP::Hep3Vector                  _n3RailCenter;
    CLHEP::Hep3Vector                  _s3RailCenter;
    // Bearing blocks on rails
    std::vector<double>                _uOutlineBBlock;
    std::vector<double>                _vOutlineBBlock;
    std::string                        _BBlockMaterial;
    double                             _lengthBBlock2;
    double                             _lengthBBlock3;
    std::vector<CLHEP::Hep3Vector>     _BBlockCenters2;
    std::vector<CLHEP::Hep3Vector>     _BBlockCenters3;

    // MBS spherical support structure
    bool                _hasMBSS;
    double              _lengthMBSS;
    std::vector<double> _uOutlineMBSS;
    std::vector<double> _vOutlineMBSS;
    CLHEP::Hep3Vector   _locationMBSS;
    std::string         _materialMBSS;

    // Cable Runs
    bool                _hasCableRunCal;
    double              _lengthCableRunCal;
    double              _rInCableRunCal;
    double              _rOutCableRunCal;
    double              _zCCableRunCal;
    double              _phi0CableRunCal;
    double              _dPhiCableRunCal;
    std::string         _materialCableRunCal;
    
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
