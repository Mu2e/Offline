//

#include "CalPatRec/inc/CalHelixPoint.hh"


namespace mu2e {

//-----------------------------------------------------------------------------
// CalHelixPoint
//-----------------------------------------------------------------------------
  // statics
  double       CalHelixPoint::_efac(1.0);
  StrawHitFlag CalHelixPoint::_useflag;

  CalHelixPoint::CalHelixPoint(size_t index, 
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

  CalHelixPoint::CalHelixPoint(size_t ind, 
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
  

  bool CalHelixPoint::use() const { 
    return !_flag.hasAnyProperty(_useflag);
  }

  bool CalHelixPoint::stereo() const {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    return _flag.hasAllProperties(stereo);
  }

  void CalHelixPoint::setUse(bool use) {
    static StrawHitFlag other(StrawHitFlag::other);
    if(!use)
      _flag.merge(other);
    else
      _flag.clear(other);
  }

  void CalHelixPoint::setOutlier(){
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    _flag.merge(outlier);
  }

  bool CalHelixPoint::isOutlier() const {
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    return _flag.hasAllProperties(outlier);
  }

  bool CalHelixPoint::isCalosel() const {
    static StrawHitFlag calosel(StrawHitFlag::calosel);
    return _flag.hasAllProperties(calosel);
  }

}
