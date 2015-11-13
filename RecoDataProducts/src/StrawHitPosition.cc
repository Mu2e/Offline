//
// Class to describe derived information from a StrawHit, in particular position.
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "TrackerGeom/inc/Tracker.hh"

namespace {
  constexpr const double invsqrt12 = 1.0/sqrt(12.0);
}

namespace mu2e {

  StrawHitPosition::StrawHitPosition(StrawHit const& hit,Straw const& straw ,SHInfo const& shinfo, StrawHitFlag const flag ) : _pos(shinfo._pos),_wdist(shinfo._tddist),
  _stindex(-1),_flag(flag) {
    CLHEP::Hep3Vector rhat = _pos.perpPart().unit();
    CLHEP::Hep3Vector phat(-rhat.y(),rhat.x(),0.0);
    double sres = 2*invsqrt12*straw.getRadius();
    double pcos = fabs(straw.getDirection().dot(phat));
    double rcos = fabs(straw.getDirection().dot(rhat));
    _pres = shinfo._tdres*pcos + sres*rcos;
    _rres = sres*pcos + shinfo._tdres*rcos;
  }

  StrawHitPosition::StrawHitPosition(StrawHitPosition const& shpos, StrawHitFlag const orflag) :
    _pos(shpos._pos),_wdist(shpos._wdist),_pres(shpos._pres),_rres(shpos._rres),_stindex(shpos._stindex),_flag(shpos._flag)
  {
    _flag.merge(orflag);
  }

  StrawHitPosition::StrawHitPosition(StereoHitCollection const& sthits, size_t stindex, size_t shindex) :
    _pos(sthits.at(stindex).pos()),
    _pres(invsqrt12*sthits.at(stindex).dist()),
    _rres(invsqrt12*sthits.at(stindex).dist()),_stindex(stindex),_flag(StrawHitFlag::stereo)
  {
    _wdist = (shindex==sthits.at(stindex).hitIndex1())? sthits.at(stindex).wdist1() : sthits.at(stindex).wdist2();
  }

  StrawHitPosition::StrawHitPosition() :_pres(-1.0),_rres(-1.0),_stindex(-1) {}

  StrawHitPosition::~StrawHitPosition() {}

  StrawHitPosition& StrawHitPosition::operator =(StrawHitPosition const& other) {
    if(this != &other){
      _pos = other._pos;
      _wdist = other._wdist;
      _pres = other._pres;
      _rres = other._rres;
      _stindex = other._stindex;
      _flag = other._flag;
    }
    return *this;
  }
}


