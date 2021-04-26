#ifndef PTMonGeom_PTMon_hh
#define PTMonGeom_PTMon_hh

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "PTMonGeom/inc/PTMonPWC.hh"

// ProductionTarget Monitor (PTMon) Object
//
// Author: Helenka Casler
//

namespace mu2e {

  class PTMon : virtual public Detector {

  public:
    PTMon(CLHEP::Hep3Vector const& originInMu2e, 
          CLHEP::HepRotation const& rotationInMu2e, 
          PTMonPWC const& nearPWC, PTMonPWC const& farPWC,
          double pwcSeparation); // TODO: better way to add/store these than just pointers
    PTMon() {}

    CLHEP::Hep3Vector const &  originInMu2e()   const { return _originInMu2e; }
    CLHEP::HepRotation const & rotationInMu2e() const { return _rotationInMu2e; }

    const PTMonPWC* nearPWC() const { return _nearPWC.get(); }
    const PTMonPWC* farPWC()  const { return _farPWC.get(); }

    double totalHeight() const { return _totalHeight; }
    double totalWidth()  const { return _totalWidth; }
    double totalLength() const { return _totalLength; }



  private:
    CLHEP::Hep3Vector _originInMu2e;
    CLHEP::HepRotation _rotationInMu2e;

    std::unique_ptr<PTMonPWC> _nearPWC; // TODO should these actually be unique_ptrs or something else?
    std::unique_ptr<PTMonPWC> _farPWC;

    double _totalHeight;
    double _totalWidth;
    double _totalLength;



  }; // class PTMon

} // namespace mu2e

#endif