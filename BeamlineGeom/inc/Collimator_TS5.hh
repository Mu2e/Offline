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
 
    double rIn()         const { return _rIn;         }
    double rMid1()       const { return _rMid1;       }
    double rMid2()       const { return _rMid2;       }
    double rOut()        const { return _rOut;        }
    double halfLengthU() const { return _halfLengthU; }
    double halfLengthD() const { return _halfLengthD; }

    std::string material()    const { return _material;    }
    std::string absMaterial() const { return _absMaterial; }

  private:
    
    double _rIn;
    double _rMid1;
    double _rMid2;
    double _rOut;
    double _halfLengthU;
    double _halfLengthD;

    std::string _material;
    std::string _absMaterial;

};

}
#endif /* BeamlineGeom_CollimatorTS5_hh */
