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

    //----
    // Primary inputs: sizes

    const std::vector<double>& coreHalfSize() const { return _coreHalfSize; }

    const std::vector<double>& neutronCaveHalfSize() const { return _neutronCaveHalfSize; }

    const std::vector<double>& mouthHalfSize() const { return _mouthHalfSize; }

    const std::vector<double>& magnetPitHalfSize() const { return _magnetPitHalfSize; }

    double minCoreShieldingThickness() const { return _minCoreShieldingThickness; }

    double collimator1horizontalLength() const { return _collimator1horizLength; }
    double collimator2horizontalLength() const { return _collimator2horizLength; }

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

    double _collimator1horizLength;
    double _collimator2horizLength;

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
