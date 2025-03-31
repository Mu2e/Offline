#ifndef BeamlineGeom_CollimatorTS1_hh
#define BeamlineGeom_CollimatorTS1_hh

//
// Class to represent the collimators. Position is relative to
// the corresponding TS straight section
//
#include "Offline/BeamlineGeom/inc/Collimator.hh"
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

  class CollimatorTS1 : public Collimator {

  friend class BeamlineMaker;

  public:

    CollimatorTS1() : Collimator() {}

    // use compiler-generated copy c'tor, copy assignment, and d'tor
    CollimatorTS1(double halfZ, CLHEP::Hep3Vector const& origin) :
      Collimator( halfZ, origin) {}

    double rIn1() const { return _rIn1; }
    double rIn2() const { return _rIn2; }
    double rIn3() const { return _rIn3; }
    double rIn4() const { return _rIn4; }
    double rOut() const { return _rOut4; }
    double rOut1()const { return _rOut1; }

    std::string material1() const { return _material1; }
    std::string material2() const { return _material2; }
    std::string material3() const { return _material3; }

    double collarHalfLength() const { return _collarHalfLength; }
    double collarZ() const { return _collarZ; }
    double collarMarginZ() const { return _collarMarginZ; }
    double collarrIn() const { return _collarrIn; }
    double collarphiBegin() const { return _collarphiBegin; }
    double collarphiDelta() const { return _collarphiDelta; }

  private:

    double _rIn1;
    double _rIn2;
    double _rIn3;
    double _rIn4;
    double _rOut4;
    double _rOut1;

    std::string _material1;
    std::string _material2;
    std::string _material3;

    double _collarHalfLength;
    double _collarZ;
    double _collarMarginZ;
    double _collarrIn;
    double _collarphiBegin;
    double _collarphiDelta;

};

}
#endif /* BeamlineGeom_CollimatorTS1_hh */
