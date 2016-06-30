#ifndef BeamlineGeom_CollimatorTS5_hh
#define BeamlineGeom_CollimatorTS5_hh

//
// Class to represent the collimators. Position is relative to
// the corresponding TS straight section
//
#include "BeamlineGeom/inc/Collimator.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class CollimatorTS5 : public Collimator {

  friend class BeamlineMaker;

  public:

    CollimatorTS5() : Collimator() {}

    // use compiler-generated copy c'tor, copy assignment, and d'tor
    CollimatorTS5(double halfZ, CLHEP::Hep3Vector origin) :  
      Collimator( halfZ, origin) {}
 
    double rIn()        const { return _rIn;        }
    double rOut()       const { return _rOut;       }

    std::string material()    const { return _material;    }

  private:
    
    double _rIn;
    double _rOut;

    std::string _material;

};

}
#endif /* BeamlineGeom_CollimatorTS5_hh */
