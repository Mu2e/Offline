#ifndef PTMGeom_PTMHead_hh
#define PTMGeom_PTMHead_hh

#include <memory>

#include "CLHEP/Vector/Rotation.h"
#include "CLHEP/Vector/ThreeVector.h"

#include "Offline/PTMGeom/inc/PTMPWC.hh"

// ProductionTarget Monitor Head (PTM) Object
//
// Author: Helenka Casler
//

// TODO: ADD THE FOLLOWING:
// - aluminum frame holding the PWCs apart
// - handle on top of the frame
// - brackets that let it sit on the top of the stand

namespace mu2e {

  class PTMHead{

  public:
    PTMHead(CLHEP::Hep3Vector const& originInMu2e, 
          CLHEP::HepRotation const& rotationInMu2e, 
          std::shared_ptr<PTMPWC> nearPWC, 
          std::shared_ptr<PTMPWC> farPWC,
          double pwcSeparation,
          std::shared_ptr<Box> holderExtrusionLong,
          std::shared_ptr<Box> holderExtrusionShort,
          std::string holderExtrusionMaterialName,
          double holderExtrusionLongSep,
          double holderExtrusionShortPos,
          double motherMargin,
          CLHEP::Hep3Vector const& parentOriginInMu2e, 
          CLHEP::HepRotation const& parentRotationInMu2e);
    PTMHead() {}

    CLHEP::Hep3Vector const &  originInMu2e()   const { return _originInMu2e; }
    CLHEP::HepRotation const & rotationInMu2e() const { return _rotationInMu2e; }
    CLHEP::Hep3Vector const &  originInParent()   const { return _originInParent; }
    CLHEP::HepRotation const & rotationInParent() const { return _rotationInParent; }

    const PTMPWC* nearPWC() const { return _nearPWC.get(); }
    const PTMPWC* farPWC()  const { return _farPWC.get(); }

    const Box* holderExtrusionLong()   const { return _holderExtrusionLong.get(); }
    const Box* holderExtrusionShort()  const { return _holderExtrusionShort.get(); }
    std::string holderExtrusionMaterialName() const { return _holderExtrusionMaterialName; }
    double holderExtrusionLongSep()    const { return _holderExtrusionLongSep; }
    double holderExtrusionShortPos()   const { return _holderExtrusionShortPos; }

    double totalHeight() const { return _totalHeight; }
    double totalWidth()  const { return _totalWidth; }
    double totalLength() const { return _totalLength; }




  private:
    CLHEP::Hep3Vector _originInMu2e;
    CLHEP::HepRotation _rotationInMu2e;
    CLHEP::Hep3Vector _originInParent;
    CLHEP::HepRotation _rotationInParent;

    std::shared_ptr<PTMPWC> _nearPWC;
    std::shared_ptr<PTMPWC> _farPWC;

    // Al extrusions holding the PWC's in position relative to each other
    std::shared_ptr<Box> _holderExtrusionLong;
    std::shared_ptr<Box> _holderExtrusionShort;
    std::string _holderExtrusionMaterialName;
    double _holderExtrusionLongSep;
    double _holderExtrusionShortPos;

    double _totalHeight;
    double _totalWidth;
    double _totalLength;



  }; // class PTMHead

} // namespace mu2e

#endif
