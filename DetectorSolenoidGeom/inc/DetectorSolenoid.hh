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

#include "Offline/Mu2eInterfaces/inc/Detector.hh"

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

    // Experimental inner lining on cryostat, requested by Jim Miller for
    // studies of shielding.
    bool hasInnerLining() const { return _hasInnerLining; }
    double innerLiningThickness() const { return _innerLiningThickness; }
    std::string innerLiningMaterial() const { return _innerLiningMaterial; }


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
    // Couplers between the bearing blocks
    double widthCoupler()  const  { return _widthCoupler; }
    double heightCoupler() const  { return _heightCoupler; }
    double yCenterCoupler() const { return _yCenterCoupler; }
    int    couplerScheme() const  { return _couplerScheme; }
    // 0 = couplers both rails, 1 = north rail only, 2 = south only

    // MBS Spherical Support outline
    bool hasMBSS() const { return _hasMBSS; }
    double  lengthMBSS()               const { return _lengthMBSS; }
    std::vector<double> uOutlineMBSS() const { return _uOutlineMBSS; }
    std::vector<double> vOutlineMBSS() const { return _vOutlineMBSS; }
    CLHEP::Hep3Vector MBSSlocation()   const { return _locationMBSS; }
    std::string MBSSmaterial()         const { return _materialMBSS; }

    // Cable Runs
    int  cableRunVersion()             const { return _cableRunVersion; }
    bool hasCableRunCal()              const { return _hasCableRunCal; }
    std::string calCableRunMaterial()  const { return _materialCableRunCal; }
    double upRInCableRunCal()          const { return _upRInCableRunCal; }
    double upROutCableRunCal()         const { return _upROutCableRunCal;}
    double upHL2CableRunCal()          const { return _upHL2CableRunCal; }
    double upZC2CableRunCal()          const { return _upZC2CableRunCal; }
    double rInCableRunCal()            const { return _rInCableRunCal; }
    double rOutCableRunCal()           const { return _rOutCableRunCal; }
    double lengthCableRunCal()         const { return _lengthCableRunCal; }
    double phi0CableRunCal()           const { return _phi0CableRunCal; }
    double dPhiCableRunCal()           const { return _dPhiCableRunCal; }
    double zCCableRunCal()             const { return _zCCableRunCal; }

    double rCableRunCalCoreFract()     const { return _rCableRunCalCoreFract; }
    double rdCableRunCalCoreFract()    const { return _rdCableRunCalCoreFract; }
    double dPhiCableRunCalCoreFract()  const { return _dPhiCableRunCalCoreFract; }
    std::string materialCableRunCalCore()   const { return _materialCableRunCalCore; }

    bool hasCableRunTrk()              const { return _hasCableRunTrk; }
    std::string trkCableRunMaterial()  const { return _materialCableRunTrk; }
    double rInCableRunTrk()            const { return _rInCableRunTrk; }
    double rOutCableRunTrk()           const { return _rOutCableRunTrk; }
    double lengthCableRunTrk()         const { return _lengthCableRunTrk; }
    double phi0CableRunTrk()           const { return _phi0CableRunTrk; }
    double dPhiCableRunTrk()           const { return _dPhiCableRunTrk; }
    double zCCableRunTrk()             const { return _zCCableRunTrk; }
    //Cabling outside IFB
    double calR1CableRunIFB()          const { return _calR1CableRunIFB   ;}
    double calR2CableRunIFB()          const { return _calR2CableRunIFB   ;}
    double calPhi0CableRunIFB()        const { return _calPhi0CableRunIFB ;}
    double calDPhiCableRunIFB()        const { return _calDPhiCableRunIFB ;}
    double calREndCableRunIFB()        const { return _calREndCableRunIFB ;}
    double calEndWCableRunIFB()        const { return _calEndWCableRunIFB ;}
    double calPhiECableRunIFB()        const { return _calPhiECableRunIFB ;}
    //Calo IFB patch panel info
    double calPR1CableRunIFB()         const { return _calPR1CableRunIFB  ;}
    double calPR2CableRunIFB()         const { return _calPR2CableRunIFB  ;}
    double calPPhi0CableRunIFB()       const { return _calPPhi0CableRunIFB;}
    double calPDPhiCableRunIFB()       const { return _calPDPhiCableRunIFB;}
    double calPZInCableRunIFB()        const { return _calPZInCableRunIFB ;}
    double calPZHLCableRunIFB()        const { return _calPZHLCableRunIFB ;}
    double calPZOutCableRunIFB()       const { return _calPZOutCableRunIFB;}
    std::string calPMatCableRunIFB()   const { return _calPMatCableRunIFB ;}
    //Calo cabling at bottom of IFB cabling
    double calBCXCableRunIFB()         const { return _calBCXCableRunIFB  ;}
    double calBLCableRunIFB()          const { return _calBLCableRunIFB   ;}

    double trkR1CableRunIFB()          const { return _trkR1CableRunIFB   ;}
    double trkR2CableRunIFB()          const { return _trkR2CableRunIFB   ;}
    double trkPhi0CableRunIFB()        const { return _trkPhi0CableRunIFB ;}
    double trkDPhiCableRunIFB()        const { return _trkDPhiCableRunIFB ;}
    double trkREndCableRunIFB()        const { return _trkREndCableRunIFB ;}
    double trkEndWCableRunIFB()        const { return _trkEndWCableRunIFB ;}
    double trkPhiECableRunIFB()        const { return _trkPhiECableRunIFB ;}
    //Tracker IFB patch panel info
    double trkPR1CableRunIFB()         const { return _trkPR1CableRunIFB  ;}
    double trkPR2CableRunIFB()         const { return _trkPR2CableRunIFB  ;}
    double trkPPhi0CableRunIFB()       const { return _trkPPhi0CableRunIFB;}
    double trkPDPhiCableRunIFB()       const { return _trkPDPhiCableRunIFB;}
    double trkPZInCableRunIFB()        const { return _trkPZInCableRunIFB ;}
    double trkPZHLCableRunIFB()        const { return _trkPZHLCableRunIFB ;}
    double trkPZOutCableRunIFB()       const { return _trkPZOutCableRunIFB;}
    std::string trkPMatCableRunIFB()   const { return _trkPMatCableRunIFB ;}
    //Tracker cabling at bottom of IFB cabling
    double trkBCXCableRunIFB()         const { return _trkBCXCableRunIFB  ;}
    double trkBLCableRunIFB()          const { return _trkBLCableRunIFB   ;}

    double zHLCableRunIFB()            const { return _zHLCableRunIFB     ;}
    std::string materialCalCableRunIFB()  const { return _materialCalCableRunIFB;}
    std::string materialTrkCableRunIFB()  const { return _materialTrkCableRunIFB;}
    double zCCableRunIFB()             const { return _zCCableRunIFB      ;}

    double rCableRunTrkCoreFract()     const { return _rCableRunTrkCoreFract; }
    double rdCableRunTrkCoreFract()    const { return _rdCableRunTrkCoreFract; }
    double dPhiCableRunTrkCoreFract()  const { return _dPhiCableRunTrkCoreFract; }
    std::string materialCableRunTrkCore()   const { return _materialCableRunTrkCore; }

    // Services pipes along bottom of DS
    bool   hasServicePipes()           const { return _hasServicePipes; }
    double servicePipeRIn()            const { return _servicePipeRIn;  }
    double servicePipeROut()           const { return _servicePipeROut; }
    double servicePipeHalfLength()     const { return _servicePipeHL;   }
    std::string servicePipeMaterial()       const { return _servicePipeMat;}
    std::string servicePipeFillMat()        const { return _servicePipeFillMat; }
    double servicePipeZC()             const { return _servicePipeZC;  }
    double servicePipeYC()             const { return _servicePipeYC;  }
    std::vector<double> servicePipeXCs() const { return _servicePipeXCs;}

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

    // experimental innerLining for shielding studies, req by Jim Miller
    bool _hasInnerLining;
    double _innerLiningThickness;
    std::string _innerLiningMaterial;

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
    // Couplers between bearing blocks
    double                             _widthCoupler;
    double                             _heightCoupler;
    double                             _yCenterCoupler;
    int                                _couplerScheme;
    // 0 = couplers along both rails,
    // 1 = couplers along north rails only,
    // 2 = couplers along south rail only

    // MBS spherical support structure
    bool                _hasMBSS;
    double              _lengthMBSS;
    std::vector<double> _uOutlineMBSS;
    std::vector<double> _vOutlineMBSS;
    CLHEP::Hep3Vector   _locationMBSS;
    std::string         _materialMBSS;

    // Cable Runs
    int                 _cableRunVersion;
    bool                _hasCableRunCal;
    double              _lengthCableRunCal;
    double              _upRInCableRunCal;
    double              _upROutCableRunCal;
    double              _upHL2CableRunCal;
    double              _upZC2CableRunCal;
    double              _rInCableRunCal;
    double              _rOutCableRunCal;
    double              _zCCableRunCal;
    double              _phi0CableRunCal;
    double              _dPhiCableRunCal;
    std::string         _materialCableRunCal;
    double              _rCableRunCalCoreFract;
    double              _rdCableRunCalCoreFract;
    double              _dPhiCableRunCalCoreFract;
    std::string         _materialCableRunCalCore;

    bool                _hasCableRunTrk;
    double              _lengthCableRunTrk;
    double              _rInCableRunTrk;
    double              _rOutCableRunTrk;
    double              _zCCableRunTrk;
    double              _phi0CableRunTrk;
    double              _dPhiCableRunTrk;
    std::string         _materialCableRunTrk;

    double              _rCableRunTrkCoreFract;
    double              _rdCableRunTrkCoreFract;
    double              _dPhiCableRunTrkCoreFract;
    std::string         _materialCableRunTrkCore;

    //Cabling outside IFB
    double              _calR1CableRunIFB   ;
    double              _calR2CableRunIFB   ;
    double              _calPhi0CableRunIFB ;
    double              _calDPhiCableRunIFB ;
    double              _calREndCableRunIFB ;
    double              _calEndWCableRunIFB ;
    double              _calPhiECableRunIFB ;
    double              _calPR1CableRunIFB  ;
    double              _calPR2CableRunIFB  ;
    double              _calPPhi0CableRunIFB;
    double              _calPDPhiCableRunIFB;
    double              _calPZInCableRunIFB ;
    double              _calPZHLCableRunIFB ;
    double              _calPZOutCableRunIFB;
    std::string         _calPMatCableRunIFB ;
    double              _calBCXCableRunIFB  ;
    double              _calBLCableRunIFB   ;

    double              _trkR1CableRunIFB   ;
    double              _trkR2CableRunIFB   ;
    double              _trkPhi0CableRunIFB ;
    double              _trkDPhiCableRunIFB ;
    double              _trkR1EndCableRunIFB;
    double              _trkREndCableRunIFB ;
    double              _trkEndWCableRunIFB ;
    double              _trkPhiECableRunIFB ;
    double              _trkPR1CableRunIFB  ;
    double              _trkPR2CableRunIFB  ;
    double              _trkPPhi0CableRunIFB;
    double              _trkPDPhiCableRunIFB;
    double              _trkPZInCableRunIFB ;
    double              _trkPZHLCableRunIFB ;
    double              _trkPZOutCableRunIFB;
    std::string         _trkPMatCableRunIFB ;
    double              _trkBCXCableRunIFB  ;
    double              _trkBLCableRunIFB   ;

    double              _zHLCableRunIFB     ;
    std::string         _materialCalCableRunIFB;
    std::string         _materialTrkCableRunIFB;
    double              _zCCableRunIFB      ;

    // Service pipes
    bool                  _hasServicePipes;
    double                _servicePipeRIn;
    double                _servicePipeROut;
    double                _servicePipeHL;
    std::string           _servicePipeMat;
    std::string           _servicePipeFillMat;
    double                _servicePipeZC;
    double                _servicePipeYC;
    std::vector<double>   _servicePipeXCs;


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
