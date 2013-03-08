//
// Class to describe derived information from a StrawHit, in particular position.
// 
// $Id: StrawHitPosition.cc,v 1.2 2013/03/08 19:13:27 brownd Exp $
// $Author: brownd $
// $Date: 2013/03/08 19:13:27 $
//
// Original author David Brown
//
// Mu2e includes
#include "RecoDataProducts/inc/StrawHitPosition.hh"
#include "TrackerGeom/inc/Tracker.hh"

namespace mu2e {
// initialize static
  StrawHitFlag StrawHitPosition::_nullflag;
  double StrawHitPosition::_invsqrt12(1.0/sqrt(12.0));

  StrawHitPosition::StrawHitPosition(StrawHit const& hit, Tracker const& tracker, ConditionsHandle<TrackerCalibrations>& tcal, StrawHitFlag const& flag ){
// get information from conditions service
    SHInfo shinfo;
    Straw const& straw = tracker.getStraw(hit.strawIndex());
    tcal->StrawHitInfo(straw,hit,shinfo);
    _pos = shinfo._pos;
    _wdist = shinfo._tddist;
    CLHEP::Hep3Vector rhat = _pos.perpPart().unit();
    CLHEP::Hep3Vector phat(-rhat.y(),rhat.x(),0.0);
    double sres = 2*_invsqrt12*straw.getRadius();
    double pcos = fabs(straw.getDirection().dot(phat)); 
    double rcos = fabs(straw.getDirection().dot(rhat)); 
    _pres = shinfo._tdres*pcos + sres*rcos;
    _rres = sres*pcos + shinfo._tdres*rcos;  
    _chisq = 0.0;
    _flag = flag;
  }

  StrawHitPosition::StrawHitPosition(StrawHitPosition const& shpos,StrawHitFlag const& orflag) :
    _pos(shpos._pos),_wdist(shpos._wdist),_pres(shpos._pres),_rres(shpos._rres),_chisq(shpos._chisq),_flag(shpos._flag)
  {
    _flag.merge(orflag);
  }

  StrawHitPosition::StrawHitPosition(StereoHit const& sthit, size_t index) :
    _pos(sthit.pos()),_pres(_invsqrt12*sthit.dist()),_rres(_invsqrt12*sthit.dist()),_chisq(sthit.chisq()),_flag(StrawHitFlagDetail::stereo)
  {
    _wdist = (index==sthit.hitIndex1())? sthit.wdist1() : sthit.wdist2();
  }

  StrawHitPosition::StrawHitPosition() :_pres(-1.0),_rres(-1.0),_chisq(-1.0) {}

  StrawHitPosition::~StrawHitPosition() {}

  StrawHitPosition& StrawHitPosition::operator =(StrawHitPosition const& other) {
    if(this != &other){
      _pos = other._pos;
      _pres = other._pres;
      _rres = other._rres;
      _chisq = other._chisq;
      _flag = other._flag;
    }
    return *this;
  }
}


