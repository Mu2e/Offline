//
//  Information about position for helix finding
//  Original author: Dave Brown (LBNL)  2014
//
#include "TrkReco/inc/XYZP.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "GeometryService/inc/getTrackerOrThrow.hh"
#include "TrackerGeom/inc/Tracker.hh"
using CLHEP::Hep3Vector;
namespace mu2e {
// statics
  double XYZP::_efac(1.0);
  StrawHitFlag XYZP::_useflag;
  StrawHitFlag XYZP::_dontuseflag;
  int XYZP::_debug(0);

  XYZP::XYZP(size_t index,StrawHit const& sh, StrawHitPosition const& shp,Straw const& straw) :
    _ind(index), _pos(shp.pos()), _phi(shp.pos().phi()), _flag(shp.flag()), _wdir(straw.getDirection()),
    _perr(_efac*shp.posRes(StrawHitPosition::phi)),_rerr(_efac*shp.posRes(StrawHitPosition::rho))
  {
    static const Hep3Vector _zdir(0.0,0.0,1.0);
    _sdir = _zdir.cross(_wdir);
  }

  XYZP::XYZP(Hep3Vector const& pos, double size) :  _ind(-1),_pos(pos),
  _phi(0.0),
  _flag(StrawHitFlag::stereo),
  _wdir(Hep3Vector(1,0,0)),_sdir(Hep3Vector(0,1,0)),
  _perr(size),_rerr(size) {}

  XYZP::XYZP(size_t ind,Hep3Vector const& pos, Hep3Vector const& wdir,
      double werr, double serr) :
    _ind(ind),_pos(pos),_phi(_pos.phi()),_wdir(wdir),_sdir(wdir.y(),-wdir.x(),0.0),
    _perr(_efac*werr),_rerr(_efac*serr){}

 void
  XYZP::rinfo(Hep3Vector const& center,VALERR& rad) const {
//    static const double onethird(1.0/3.0);
//    static const double invsqrt12(1./sqrt(12.0));
// average the 1-sigma radii to account for non-linear errors
    double rvec = Hep3Vector(_pos - center).perp();
//    rad._val = onethird*(rvec+rvec1+rvec2);
    rad._val = rvec;
    rad._err = _rerr;
    if(_debug > 1)std::cout << "rinfo : r = " << rad._val << " rerr = " << rad._err  << std::endl;
  }

  void
  XYZP::finfo(Hep3Vector const& center,VALERR& phi) const {
//    static const double onethird(1.0/3.0);
//    static const double invsqrt12(1./sqrt(12.0));
// average the 1-sigma radii to account for non-linear errors
    double phi0 = Hep3Vector(_pos - center).phi();
//    rad._val = onethird*(rvec+rvec1+rvec2);
    phi._val = phi0;
    phi._err = _perr; 
    if(_debug > 1)std::cout << "finfo : phi = " << phi._val << " ferr = " << phi._err << std::endl;
 }

  bool 
  XYZP::use() const {
    return (!_flag.hasAnyProperty(_dontuseflag))
      && (_flag.hasAllProperties(_useflag) || _useflag.empty());
  }

  bool 
  XYZP::stereo() const {
    static StrawHitFlag stereo(StrawHitFlag::stereo);
    return _flag.hasAllProperties(stereo);
  }

  void XYZP::setUse(bool use) {
    static StrawHitFlag other(StrawHitFlag::other);
    if(!use)
      _flag.merge(other);
    else
      _flag.clear(other);
  }

  void
  XYZP::setOutlier(){
    static StrawHitFlag outlier(StrawHitFlag::outlier);
    _flag.merge(outlier);
  }

  void XYZP::fillXYZP(StrawHitCollection const& shcol,
    StrawHitPositionCollection const& shpcol, std::vector<size_t> hits, XYZPVector& xyzp) {
    const Tracker& tracker = getTrackerOrThrow();
    // loop over straw hits, and store their positions
    for(auto const& istr : hits) { 
      StrawHit const& sh = shcol.at(istr);
      Straw const& straw= tracker.getStraw(sh.strawIndex());
      StrawHitPosition const& shp = shpcol.at(istr);
      XYZP pos(istr,sh,shp,straw);
      xyzp.push_back(pos);
    } 
  }

}

