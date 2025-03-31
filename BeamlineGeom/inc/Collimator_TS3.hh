#ifndef BeamlineGeom_CollimatorTS3_hh
#define BeamlineGeom_CollimatorTS3_hh

//
// Class to represent the collimators. Position is relative to
// the corresponding TS straight section
//
#include "Offline/BeamlineGeom/inc/Collimator.hh"
#include "CLHEP/Vector/ThreeVector.h"

namespace mu2e {

  class CollimatorTS3 : public Collimator {

  friend class BeamlineMaker;

  public:

    CollimatorTS3() : Collimator() {}

    // use compiler-generated copy c'tor, copy assignment, and d'tor
    CollimatorTS3(double halfZ, CLHEP::Hep3Vector const& origin) :
      Collimator( halfZ, origin) {}

    double      rOut() const { return _rOut; }

    double      hole()             const { return _hole;             }
    double      holeRadius()       const { return _holeRadius;       }
    double      holeHalfHeight()   const { return _holeHalfHeight;   }
    double      holeDisplacement() const { return _holeDisplacement; }
    double      rotationAngle()    const { return _rotationAngle;    }

    bool        useFlashBlock()      const { return _useFlashBlock;    }
    double      flashBlockHeight()   const { return _flashBlockHeight; }
    double      flashBlockWidth()    const { return _flashBlockWidth;  }
    double      flashBlockLength()   const { return _flashBlockLength; }
    double      flashBlockTranOff()  const { return _flashBlockTO;     }
    double      flashBlockLongOff()  const { return _flashBlockLO;     }
    std::string flashBlockMaterial() const { return _flashBlockMaterial;}

    std::string material()           const { return _material; }

  private:

    double _rOut;

    double _hole;
    double _holeRadius;
    double _holeHalfHeight;
    double _holeDisplacement;
    double _rotationAngle;

    bool   _useFlashBlock;
    double _flashBlockHeight;
    double _flashBlockWidth;
    double _flashBlockLength;
    double _flashBlockTO;
    double _flashBlockLO;
    std::string _flashBlockMaterial;

    std::string _material;

};

}
#endif /* BeamlineGeom_CollimatorTS3_hh */
