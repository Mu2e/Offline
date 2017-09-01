//

#include "CalPatRec/inc/XYZPHack.hh"


namespace mu2e {

//-----------------------------------------------------------------------------
// XYZPHack
//-----------------------------------------------------------------------------
  // statics
  double       XYZPHack::_efac(1.0);
  StrawHitFlag XYZPHack::_useflag;

  XYZPHack::XYZPHack(size_t index, 
		     StrawHit const& sh, 
		     StrawHitPosition const& shp, 
		     Straw const& straw,
		     StrawHitFlag const& flag) :
    _ind (index), 
    _pos (shp.pos()), 
    _phi (shp.pos().phi()), 
    _flag(flag),
    _used(0),
    _wdir(straw.getDirection()),
    _straw(&straw),
    _strawhit(&sh),
    _perr(_efac*shp.posRes(StrawHitPosition::wire)),
    _rerr(_efac*shp.posRes(StrawHitPosition::trans))
  {
    static const CLHEP::Hep3Vector _zdir(0.0,0.0,1.0);
    _sdir = _zdir.cross(_wdir);
  }

  XYZPHack::XYZPHack(size_t ind, 
		     CLHEP::Hep3Vector const& pos, 
		     CLHEP::Hep3Vector const& wdir, 
		     double werr, 
		     double serr) :
    _ind(ind),
    _pos(pos),
    _phi(_pos.phi()),
    _flag(),
    _used(0),
    _wdir(wdir),
    _sdir(wdir.y(),-wdir.x(),0.0),
    _straw(0),
    _strawhit(0),
    _perr(_efac*werr),
    _rerr(_efac*serr)
  {
  }
  
//   void XYZPHack::rinfo(CLHEP::Hep3Vector const& center, VALERR& rad) const {
//     //    static const double onethird(1.0/3.0);
//     //    static const double invsqrt12(1./sqrt(12.0));
//     // average the 1-sigma radii to account for non-linear errors
//     double rvec = CLHEP::Hep3Vector(_pos - center).perp();
//     //    rad._val = onethird*(rvec+rvec1+rvec2);
//     rad._val = rvec;
//     rad._err = _rerr;
    
//   }

//   void XYZPHack::finfo(CLHEP::Hep3Vector const& center,VALERR& phi) const {
//     //    static const double onethird(1.0/3.0);
//     //    static const double invsqrt12(1./sqrt(12.0));
//     // average the 1-sigma radii to account for non-linear errors
//     double phi0 = CLHEP::Hep3Vector(_pos - center).phi();
//     //    rad._val = onethird*(rvec+rvec1+rvec2);
//     phi._val = phi0;
//     phi._err = _perr; 
//   }

  bool XYZPHack::use() const { 
    return !_flag.hasAnyProperty(_useflag);
  }

  bool XYZPHack::stereo() const {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    return _flag.hasAllProperties(stereo);
  }

  void XYZPHack::setUse(bool use) {
    static StrawHitFlag other(StrawHitFlag::other);
    if(!use)
      _flag.merge(other);
    else
      _flag.clear(other);
  }

  void XYZPHack::setOutlier(){
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    _flag.merge(outlier);
  }

  bool XYZPHack::isOutlier() const {
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    return _flag.hasAllProperties(outlier);
  }

  bool XYZPHack::isCalosel() const {
    static StrawHitFlag calosel(StrawHitFlag::calosel);
    return _flag.hasAllProperties(calosel);
  }

}
