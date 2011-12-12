// Geometry of the proton beam dump.
// 
// Andrei Gaponenko, 2011

#ifndef PROTONBEAMDUMP_HH
#define PROTONBEAMDUMP_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "GeometryService/inc/Detector.hh"

namespace mu2e {
    
  class ProtonBeamDumpMaker;

  class ProtonBeamDump : public Detector {
  public: 
    // implement Detector's method
    virtual std::string name() const { return "ProtonBeamDump"; }

    class CollimatorExtMonFNAL {
      friend class ProtonBeamDumpMaker;

      std::string _name;
      double _horizontalLength;
      std::vector<double> _channelWidth;
      std::vector<double> _channelHeight;
      std::vector<double> _alignmentPlugRadius;
      std::vector<double> _alignmentHoleRClearance;
      double _radiusTransitiondZ;
      double _angleH;
      double _angleV;
      double _entranceOffsetX;
      double _entranceOffsetY;

    public:
      std::string name() const { return _name; }
      double horizontalLength() const { return _horizontalLength; }
      const std::vector<double> &channelWidth() const { return _channelWidth; }
      const std::vector<double> &channelHeight() const { return _channelHeight; }
      const std::vector<double> &alignmentPlugRadius() const { return _alignmentPlugRadius; }
      const std::vector<double> &alignmentHoleRClearance() const { return _alignmentHoleRClearance; }
      double radiusTransitiondZ() const { return _radiusTransitiondZ; }
      double angleH() const { return _angleH; }
      double angleV() const { return _angleV; }

      // These are coordinates in the beam dump face plane
      double entranceOffsetX() const { return _entranceOffsetX; }
      double entranceOffsetY() const { return _entranceOffsetY; }
    };

    //----
    // Primary inputs: sizes

    const std::vector<double>& coreHalfSize() const { return _coreHalfSize; }

    const std::vector<double>& neutronCaveHalfSize() const { return _neutronCaveHalfSize; }

    const std::vector<double>& mouthHalfSize() const { return _mouthHalfSize; }

    const std::vector<double>& magnetPitHalfSize() const { return _magnetPitHalfSize; }

    double minCoreShieldingThickness() const { return _minCoreShieldingThickness; }

    const CollimatorExtMonFNAL& collimator1() const { return _collimator1; }
    const CollimatorExtMonFNAL& collimator2() const { return _collimator2; }

    //----
    // Primary inputs: placement

    const CLHEP::Hep3Vector& coreCenterInMu2e() const { return _coreCenterInMu2e; }
    // absolute w.r.t to the Mu2e
    double coreRotY() const { return _coreRotY; }

    //---
    // Derived stuff

    // The box containing the dump, the filter magnet room, and the two collimators 
    // for the secondaries.
    // Does not include the wedge. 
    const std::vector<double>& enclosureHalfSize() const { return _enclosureHalfSize; }
    CLHEP::Hep3Vector const& enclosureCenterInMu2e() const { return _enclosureCenterInMu2e; }
    CLHEP::HepRotation const& enclosureRotationInMu2e() const { return _enclosureRotationInMu2e; }

    CLHEP::Hep3Vector const& coreCenterInEnclosure() const { return _coreCenterInEnclosure; }
    CLHEP::Hep3Vector const& magnetPitCenterInEnclosure() const { return _magnetPitCenterInEnclosure; }

    double shieldingFaceXmin() const { return _shieldingFaceXmin; }
    double shieldingFaceXmax() const { return _shieldingFaceXmax; }

    double shieldingFaceZatXmin() const { return _shieldingFaceZatXmin; }
    double shieldingFaceZatXmax() const { return _shieldingFaceZatXmax; }

    //----------------------------------------------------------------
  private: 
    friend class ProtonBeamDumpMaker;

    // Private ctr: the class should be only obtained via ProtonBeamDumpFNAL::ProtonBeamDumpMaker.
    ProtonBeamDump() : _enclosureRotationInMu2e(CLHEP::HepRotation::IDENTITY) {}

    std::vector<double> _coreHalfSize;
    std::vector<double> _neutronCaveHalfSize;
    std::vector<double> _mouthHalfSize;
    std::vector<double> _magnetPitHalfSize;
    double _minCoreShieldingThickness;

    CollimatorExtMonFNAL _collimator1;
    CollimatorExtMonFNAL _collimator2;

    CLHEP::Hep3Vector _coreCenterInMu2e;
    double _coreRotY;
    
    // computed stuff
    std::vector<double> _enclosureHalfSize;
    CLHEP::Hep3Vector _enclosureCenterInMu2e;
    CLHEP::HepRotation _enclosureRotationInMu2e;

    CLHEP::Hep3Vector _coreCenterInEnclosure;
    CLHEP::Hep3Vector _magnetPitCenterInEnclosure;
    
    double _shieldingFaceXmin;
    double _shieldingFaceXmax;
    double _shieldingFaceZatXmin;
    double _shieldingFaceZatXmax;
  };
}

#endif/*PROTONBEAMDUMP_HH*/
