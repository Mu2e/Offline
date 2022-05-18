#ifndef PTMGeom_PTM_hh
#define PTMGeom_PTM_hh

#include <memory>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/PTMGeom/inc/PTMPWC.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"

// ProductionTarget Monitor (PTM) Object
//
// Author: Helenka Casler
//

namespace mu2e {

  class PTM : virtual public Detector {

  public:
    PTM(CLHEP::Hep3Vector const& originInMu2e,
          CLHEP::HepRotation const& rotationInMu2e,
          std::shared_ptr<PTMPWC> nearPWC,
          std::shared_ptr<PTMPWC> farPWC,
          double pwcSeparation,
          double motherMargin);
    PTM() {}

    CLHEP::Hep3Vector const &  originInMu2e()   const { return _originInMu2e; }
    CLHEP::HepRotation const & rotationInMu2e() const { return _rotationInMu2e; }

    const PTMPWC* nearPWC() const { return _nearPWC.get(); }
    const PTMPWC* farPWC()  const { return _farPWC.get(); }

    double totalHeight() const { return _totalHeight; }
    double totalWidth()  const { return _totalWidth; }
    double totalLength() const { return _totalLength; }




  private:
    CLHEP::Hep3Vector _originInMu2e;
    CLHEP::HepRotation _rotationInMu2e;

    std::shared_ptr<PTMPWC> _nearPWC;
    std::shared_ptr<PTMPWC> _farPWC;

    double _totalHeight;
    double _totalWidth;
    double _totalLength;



  }; // class PTM

} // namespace mu2e

#endif
