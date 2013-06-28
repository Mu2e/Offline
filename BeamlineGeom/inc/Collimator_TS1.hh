#ifndef BeamlineGeom_CollimatorTS1_hh
#define BeamlineGeom_CollimatorTS1_hh

//
// Class to represent the collimators. Position is relative to
// the corresponding TS straight section
//
#include "BeamlineGeom/inc/Collimator.hh"
#include "CLHEP/Vector/ThreeVector.h"


namespace mu2e {

  class CollimatorTS1 : public Collimator {

  friend class BeamlineMaker;

  public:

    CollimatorTS1() : Collimator() {}

    // use compiler-generated copy c'tor, copy assignment, and d'tor
    CollimatorTS1(double halfZ, CLHEP::Hep3Vector origin) :  
      Collimator( halfZ, origin) {}

    double rIn1() const { return _rIn1; }
    double rIn2() const { return _rIn2; }
    double rIn3() const { return _rIn3; }

    std::string material1() const { return _material1; }
    std::string material2() const { return _material2; }

  private:

    double _rIn1;
    double _rIn2;
    double _rIn3;

    std::string _material1;
    std::string _material2;
};

}
#endif /* BeamlineGeom_CollimatorTS1_hh */
