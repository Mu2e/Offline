// Geometry of the TSdA
//
// Kyle Knoepfel, 2013

#ifndef BEAMLINEGEOM_TSDA_HH
#define BEAMLINEGEOM_TSDA_HH

#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class TSdAMaker;

  class TSdA : virtual public Detector {
  public:

    double r4() const { return _r4; }

    double halfLength4() const { return _halfLength4; }

    const CLHEP::Hep3Vector& position() const { return _position; }

    std::string material4() const { return _mat4; }

    int  version() const { return _version; }
    int  build  () const { return _build  ; }

    //----------------------------------------------------------------
  private:
    friend class TSdAMaker;

    // Private ctr: the class should be only obtained via TSdA::TSdAMaker.
    TSdA();

    double _r4;

    double _halfLength4;

    std::string _mat4;

    CLHEP::Hep3Vector _position;

    int  _version;
    int  _build;

    // Needed for persistency
    //    template<class T> friend class art::Wrapper;
    //    TSdA() {}
  };
}

#endif/*BEAMLINEGEOM_TSDA_HH*/
