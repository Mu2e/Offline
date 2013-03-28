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

    // Assume fixed outer radius for entire int. neutron abs.
    double rOut() const { return _rOut; }

    // ABS1 cylinder parameters
    // - assumes multiple concentric cylinders
    // - only supports concentric cylinders of constant length
    const std::vector<double>& rInAbs1Vec() const { return _rInAbs1Vec; }
    double halfLengthAbs1() const { return _halfLengthAbs1; }
    const std::vector<std::string>& materialAbs1Vec() const { return _matAbs1Vec; }

    // ABS2 cylinder parameters
    double rInAbs2()  const { return _rInAbs2; }
    std::string materialAbs2()   const { return _matAbs2; }

    // in mu2e coordinates
    const CLHEP::Hep3Vector& positionAbs1() const { return _positionAbs1; }

    //----------------------------------------------------------------
  private:
    friend class InternalNeutronAbsorberMaker;

    // Private ctr: the class should be only obtained via InternalNeutronAbsorber::InternalNeutronAbsorberMaker.
    InternalNeutronAbsorber();

    double _rOut;

    // Abs1 private members
    std::vector <double> _rInAbs1Vec;
    double _halfLengthAbs1;
    CLHEP::Hep3Vector _positionAbs1;
    std::vector<std::string> _matAbs1Vec;

    // Abs2 private members (the position of abs2 can be determined by
    // positionAbs1 + halfLengthAbs1 + halfLengthAbs2
    // - halfLengthAbs2 is constrained based on the z-extent of the DS
    //   (so it is not included as a member here)
    double _rInAbs2;
    std::string _matAbs2;

    // Needed for persistency
    //    template<class T> friend class art::Wrapper;
    //    InternalNeutronAbsorber() {}
  };
}

#endif/*INTERNALNEUTRONABSORBER_HH*/
