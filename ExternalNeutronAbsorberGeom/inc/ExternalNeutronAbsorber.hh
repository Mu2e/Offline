// Geometry of the external neutron absorber.
//
// Kyle Knoepfel, 2013

#ifndef EXTERNALNEUTRONABSORBER_HH
#define EXTERNALNEUTRONABSORBER_HH

#include <vector>

// #include "art/Persistency/Common/Wrapper.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "Mu2eInterfaces/inc/Detector.hh"

namespace mu2e {

  class ExternalNeutronAbsorberMaker;

  class ExternalNeutronAbsorber : virtual public Detector {
  public:

    double halfLengthZ()  const  { return _halfLengthZ;   }
    double halfLengthXY() const  { return _halfLengthXY;  }
    double halfThickness() const { return _halfThickness; }

    const CLHEP::Hep3Vector& position() const { return _position; }

    std::string material() const { return _material; }

    //----------------------------------------------------------------
  private:
    friend class ExternalNeutronAbsorberMaker;

    // Private ctr: the class should be only obtained via ExternalNeutronAbsorber::ExternalNeutronAbsorberMaker.
    ExternalNeutronAbsorber();

    double _halfLengthZ;
    double _halfLengthXY;
    double _halfThickness;

    std::string _material;

    CLHEP::Hep3Vector _position;

    // Needed for persistency
    //    template<class T> friend class art::Wrapper;
    //    ExternalNeutronAbsorber() {}
  };
}

#endif/*EXTERNALNEUTRONABSORBER_HH*/
