// Geometry of the proton beam dump.
//
// Andrei Gaponenko, 2011

#ifndef PROTONBEAMDUMP_HH
#define PROTONBEAMDUMP_HH

#include <vector>

#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class ProtonBeamDumpMaker;

  class ProtonBeamDump : virtual public Detector {
  public:

    //----------------------------------------------------------------
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

    public:
      std::string name() const { return _name; }
      double horizontalLength() const { return _horizontalLength; }
      const std::vector<double> &channelWidth() const { return _channelWidth; }
      const std::vector<double> &channelHeight() const { return _channelHeight; }
      const std::vector<double> &alignmentPlugRadius() const { return _alignmentPlugRadius; }
      const std::vector<double> &alignmentHoleRClearance() const { return _alignmentHoleRClearance; }
      double radiusTransitiondZ() const { return _radiusTransitiondZ; }

      // these two are not just "positioning" parameters but also affect the shape
      // thus they belong to this class
      double angleH() const { return _angleH; }
      double angleV() const { return _angleV; }
    };

    //----------------------------------------------------------------
    class FilterMagnetExtMonFNAL {
      friend class ProtonBeamDumpMaker;
      std::vector<double> _outerHalfSize;
      double _apertureWidth;
      double _apertureHeight;
      double _fieldStrength;
    public:
      const std::vector<double> &outerHalfSize() const { return _outerHalfSize; }
      double apertureWidth() const { return _apertureWidth; }
      double apertureHeight() const { return _apertureHeight; }
      double fieldStrength() const { return _fieldStrength; }

      // derived:
      double trackBendHalfAngle(double momentum) const;
    };

    //----------------------------------------------------------------
    // Primary inputs: sizes

    const std::vector<double>& coreHalfSize() const { return _coreHalfSize; }

    const std::vector<double>& neutronCaveHalfSize() const { return _neutronCaveHalfSize; }

    const std::vector<double>& mouthHalfSize() const { return _mouthHalfSize; }

    const std::vector<double>& magnetPitHalfSize() const { return _magnetPitHalfSize; }

    double minCoreShieldingThickness() const { return _minCoreShieldingThickness; }

    //----
    // Primary inputs: placement

    const CLHEP::Hep3Vector& coreCenterInMu2e() const { return _coreCenterInMu2e; }
    // absolute w.r.t to the Mu2e
    double coreRotY() const { return _coreRotY; }

    // The offsets are w.r.t. the dump core center, in the plane of the
    // dump shielding face.
    double filterEntranceOffsetX() const { return _filterEntranceOffsetX; }
    double filterEntranceOffsetY() const { return _filterEntranceOffsetY; }

    double filterAngleH() const { return _filterAngleH; }
    double filterEntranceAngleV() const { return _filterEntranceAngleV; }

    double extMonFilter_nominalMomentum() const { return _extMonFilter_nominalMomentum; }

    const FilterMagnetExtMonFNAL& filterMagnet() const { return _filterMagnet; }
    const CollimatorExtMonFNAL& collimator1() const { return _collimator1; }
    const CollimatorExtMonFNAL& collimator2() const { return _collimator2; }

    //---
    // Derived stuff

    // The box containing the dump, the filter magnet room, and the two collimators
    // for the secondaries.
    // Does not include the wedge.
    const std::vector<double>& enclosureHalfSize() const { return _enclosureHalfSize; }
    const CLHEP::Hep3Vector& enclosureCenterInMu2e() const { return _enclosureCenterInMu2e; }
    const CLHEP::HepRotation& enclosureRotationInMu2e() const { return _enclosureRotationInMu2e; }

    const CLHEP::Hep3Vector& coreCenterInEnclosure() const { return _coreCenterInEnclosure; }
    const CLHEP::Hep3Vector& magnetPitCenterInEnclosure() const { return _magnetPitCenterInEnclosure; }

    double shieldingFaceXmin() const { return _shieldingFaceXmin; }
    double shieldingFaceXmax() const { return _shieldingFaceXmax; }

    double shieldingFaceZatXmin() const { return _shieldingFaceZatXmin; }
    double shieldingFaceZatXmax() const { return _shieldingFaceZatXmax; }

    const CLHEP::Hep3Vector& collimator1CenterInEnclosure() const { return _collimator1CenterInEnclosure; }

    const CLHEP::Hep3Vector& filterMagnetCenterInEnclosure() const { return _filterMagnetCenterInEnclosure; }
    double filterMagnetAngleV() const { return _filterMagnetAngleV; }

    const CLHEP::Hep3Vector& collimator2CenterInEnclosure() const { return _collimator2CenterInEnclosure; }

    // convenience accessor for event generators
    // returns the point of intersection of collimator1 axis with the dump shielding face
    CLHEP::Hep3Vector filterEntranceInMu2e() const;

    //----------------------------------------------------------------
    // Transform to the "beam dump" coordinate system, which is centered
    // at the core center, and is rotated around the Y axis w.r.t the mu2e system

    CLHEP::Hep3Vector mu2eToBeamDump_position(const CLHEP::Hep3Vector& mu2epos) const;
    CLHEP::Hep3Vector mu2eToBeamDump_momentum(const CLHEP::Hep3Vector& mu2emom) const;

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

    CLHEP::Hep3Vector _coreCenterInMu2e;
    double _coreRotY;

    double _filterEntranceOffsetX;
    double _filterEntranceOffsetY;
    double _filterAngleH;
    double _filterEntranceAngleV;

    double _extMonFilter_nominalMomentum;

    FilterMagnetExtMonFNAL _filterMagnet;
    CollimatorExtMonFNAL _collimator1;
    CollimatorExtMonFNAL _collimator2;

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

    CLHEP::Hep3Vector _collimator1CenterInEnclosure;

    CLHEP::Hep3Vector _filterMagnetCenterInEnclosure;
    double _filterMagnetAngleV;

    CLHEP::Hep3Vector _collimator2CenterInEnclosure;
  };
}

#endif/*PROTONBEAMDUMP_HH*/
