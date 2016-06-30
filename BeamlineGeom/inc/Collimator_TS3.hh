#ifndef BeamlineGeom_CollimatorTS3_hh
#define BeamlineGeom_CollimatorTS3_hh

//
// Class to represent the collimators. Position is relative to
// the corresponding TS straight section
//
#include "BeamlineGeom/inc/Collimator.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class CollimatorTS3 : public Collimator {

  friend class BeamlineMaker;

  public:

    CollimatorTS3() : Collimator() {}

    // use compiler-generated copy c'tor, copy assignment, and d'tor
    CollimatorTS3(double halfZ, CLHEP::Hep3Vector origin) :  
      Collimator( halfZ, origin) {}

    double hole()             const { return _hole;             }
    double holeRadius()       const { return _holeRadius;       }
    double holeHalfHeight()   const { return _holeHalfHeight;   }
    double holeDisplacement() const { return _holeDisplacement; }
    double rotationAngle()    const { return _rotationAngle;    }

    std::string material()    const { return _material; }

  private:

    double _hole;
    double _holeRadius;
    double _holeHalfHeight;
    double _holeDisplacement;
    double _rotationAngle;

    std::string _material;

};

}
#endif /* BeamlineGeom_CollimatorTS3_hh */
