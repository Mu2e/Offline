// Geometry of the detector solenoid
//
// This interface accepts the following inputs:
//  - the inner and outer radii of the DS
//  - the half-length of the DS
//  - the half-lengths of the three vacuum volumes 
//    located inside of the DS volume
//  - the half-length of the front face of the solenoid (Al)
//  - the materials of the DS (default Al), and the vacuum 
//    vols (default DSVacuum)
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

    // solenoid cylinder parameters
    double rIn()  const { return _rIn; }
    double rOut()  const { return _rOut; }
    double halfLength() const { return _halfLength; }

    // in mu2e coordinates 
    const CLHEP::Hep3Vector& position() const { return _position; }

    // The subdivision of the DS vacuum volume is not physical,
    // but it needs to be in Geometry because real physical
    // pieces are placed inside.
    double halfLengthDs1() const { return _ds1HalfLength; }
    double halfLengthDs2() const { return _ds2HalfLength; }

    // The half-thickness of the front wall of the DS
    double frontHalfLength() const { return _frontHalfLength; }

    std::string material() const { return _materialName; }
    std::string insideMaterial() const { return _insideMaterialName; }


    //----------------------------------------------------------------
  private:
    friend class DetectorSolenoidMaker;

    // Private ctr: the class should be only obtained via DetectorSolenoid::DetectorSolenoidMaker.
    DetectorSolenoid();

    // solenoid features
    double _rIn;
    double _rOut;
    double _halfLength;
    CLHEP::Hep3Vector _position;
    double _ds1HalfLength;
    double _ds2HalfLength;
    double _ds3HalfLength;
    double _frontHalfLength;
    std::string _materialName;
    std::string _insideMaterialName;

    // Needed for persistency
    //    template<class T> friend class art::Wrapper;
    //    DetectorSolenoid() {}
  };
}

#endif/*DETECTORSOLENOID_HH*/
