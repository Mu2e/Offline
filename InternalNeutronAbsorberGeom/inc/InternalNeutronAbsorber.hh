// Geometry of the internal neutron absorber.
//
// Kyle Knoepfel, 2013

#ifndef INTERNALNEUTRONABSORBER_HH
#define INTERNALNEUTRONABSORBER_HH

#include <vector>

// #include "art/Persistency/Common/Wrapper.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class InternalNeutronAbsorberMaker;

  class InternalNeutronAbsorber : virtual public Detector {
  public:

    double rOut() const { return _rOut; }
    double r4() const { return _r4; }
    double rIn3() const { return _rIn3; }
    double rIn2() const { return _rIn2; }
    double rIn1() const { return _rIn1; }

    double halfLength4() const { return _halfLength4; }
    double halfLength3() const { return _halfLength3; }
    double halfLength2() const { return _halfLength2; }
    double halfLength1() const { return _halfLength1; }

    const CLHEP::Hep3Vector& position() const { return _position; }

    std::string material1() const { return _mat1; }
    std::string material2() const { return _mat2; }
    std::string material3() const { return _mat3; }
    std::string material4() const { return _mat4; }

    //----------------------------------------------------------------
  private:
    friend class InternalNeutronAbsorberMaker;

    // Private ctr: the class should be only obtained via InternalNeutronAbsorber::InternalNeutronAbsorberMaker.
    InternalNeutronAbsorber();

    double _rOut;
    double _r4;
    double _rIn3;
    double _rIn2;
    double _rIn1;

    double _halfLength4;
    double _halfLength3;
    double _halfLength2;
    double _halfLength1;

    
    std::string _mat4;
    std::string _mat3;
    std::string _mat2;
    std::string _mat1;

    CLHEP::Hep3Vector _position;

    // Needed for persistency
    //    template<class T> friend class art::Wrapper;
    //    InternalNeutronAbsorber() {}
  };
}

#endif/*INTERNALNEUTRONABSORBER_HH*/
