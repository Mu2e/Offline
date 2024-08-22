#ifndef BeamlineGeom_CollimatorTS5_hh
#define BeamlineGeom_CollimatorTS5_hh

//
// Class to represent the collimators. Position is relative to
// the corresponding TS straight section
//
#include "Offline/BeamlineGeom/inc/Collimator.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class CollimatorTS5 : public Collimator {

  friend class BeamlineMaker;

  public:

    CollimatorTS5() : Collimator() {}

    // use compiler-generated copy c'tor, copy assignment, and d'tor
    CollimatorTS5(double halfZ, CLHEP::Hep3Vector const& origin) :
      Collimator( halfZ, origin) {}

    double rIn()        const { return _rIn;        }
    double rOut()       const { return _rOut;       }
    int    version()    const { return _version;    }

    std::string material()    const { return _material;    }

  private:

    double _rIn;
    double _rOut;
    int    _version;

    std::string _material;

};

}
#endif /* BeamlineGeom_CollimatorTS5_hh */
