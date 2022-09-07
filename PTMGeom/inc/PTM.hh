#ifndef PTMGeom_PTM_hh
#define PTMGeom_PTM_hh

#include <memory>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/PTMGeom/inc/PTMPWC.hh"
#include "Offline/PTMGeom/inc/PTMStand.hh"
#include "Offline/PTMGeom/inc/PTMHead.hh"
#include "Offline/Mu2eInterfaces/inc/Detector.hh"

// ProductionTarget Monitor (PTM) Object
//
// Author: Helenka Casler
//

namespace mu2e {

  class PTM : virtual public Detector {

  public:
    PTM(int version,
          CLHEP::Hep3Vector const& originInMu2e,
          CLHEP::HepRotation const& rotationInMu2e,
          std::shared_ptr<PTMStand> ptmStand,
          std::shared_ptr<PTMHead> ptmHead);
    PTM(int version,
          CLHEP::Hep3Vector const& originInMu2e,
          CLHEP::HepRotation const& rotationInMu2e,
          std::shared_ptr<PTMPWC> nearPWC,
          std::shared_ptr<PTMPWC> farPWC,
          double pwcSeparation,
          double motherMargin);
    PTM() {}

    int version() const { return _version; }

    CLHEP::Hep3Vector const &  originInMu2e()   const { return _originInMu2e; }
    CLHEP::HepRotation const & rotationInMu2e() const { return _rotationInMu2e; }

    const PTMStand* ptmStand() const { return _ptmStand.get(); }
    const PTMHead* ptmHead()  const { return _ptmHead.get(); }

    const PTMPWC* nearPWC() const { return _nearPWC.get(); }
    const PTMPWC* farPWC()  const { return _farPWC.get(); }

    double totalHeight() const { return _totalHeight; }
    double totalWidth()  const { return _totalWidth; }
    double totalLength() const { return _totalLength; }




  private:
    int _version;

    CLHEP::Hep3Vector _originInMu2e;
    CLHEP::HepRotation _rotationInMu2e;

    std::shared_ptr<PTMStand> _ptmStand;
    std::shared_ptr<PTMHead> _ptmHead;

    std::shared_ptr<PTMPWC> _nearPWC;
    std::shared_ptr<PTMPWC> _farPWC;

    double _totalHeight;
    double _totalLength; // Z dimension in its own coordinates
    double _totalWidth;



  }; // class PTM

} // namespace mu2e

#endif
