//
// BaBar includes
#include "BTrk/BaBar/BaBar.hh"
#include "Mu2eBTrk/inc/BaBarMu2eField.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "BFieldGeom/inc/BFieldManager.hh"
#include "GeometryService/inc/DetectorSystem.hh"

namespace mu2e
{

  BaBarMu2eField::BaBarMu2eField(CLHEP::Hep3Vector const& origin) : _bnom(0.0), _origin(origin) {
  }

  BaBarMu2eField::~BaBarMu2eField(){}

  // BaBar interface.  Note we have to change units here to the BaBar conventions
  CLHEP::Hep3Vector
  BaBarMu2eField::bFieldVect (const HepPoint &point)const {
    static GeomHandle<BFieldManager> bfmgr;
    static GeomHandle<DetectorSystem> det;

    // change coordinates to mu2e
    CLHEP::Hep3Vector vpoint(point.x(),point.y(),point.z());
    CLHEP::Hep3Vector vpoint_mu2e = det->toMu2e(vpoint);
    CLHEP::Hep3Vector field = bfmgr->getBField(vpoint_mu2e);
    return field;
  }

  double
  BaBarMu2eField::bFieldNominal() const {
    if(_bnom == 0.0){
      // field at tracker origin
      // nominal field is the z component of the field at the origin
      HepPoint po(_origin.x(),_origin.y(),_origin.z());
      CLHEP::Hep3Vector bo = bFieldVect(po);
      _bnom = bo.z();
    }
    return _bnom;
  }
}


