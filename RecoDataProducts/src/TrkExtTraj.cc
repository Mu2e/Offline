//
//  $Id: TrkExtTraj.cc,v 1.1 2012/08/04 00:16:02 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/08/04 00:16:02 $
//
//  Original author MyeongJae Lee
//
//


#include "RecoDataProducts/inc/TrkExtTraj.hh"
#include "RecoDataProducts/inc/TrkExtTrajPoint.hh"

using namespace std;

namespace mu2e {

  TrkExtTraj::TrkExtTraj() {
    _pt.clear();
    _pahitidx.clear();
    _sthitidx.clear();
    _ret = -1;
  }

  TrkExtTraj::TrkExtTraj(const TrkExtTrajPoint & trajPoint) {
    _pt.clear();
    _pt.push_back(trajPoint);
    _pahitidx.clear();
    _sthitidx.clear();
    _ret = -1;
  }

  TrkExtTraj::TrkExtTraj(const TrkExtTraj & dt) {
    _pt = dt._pt;
    _ret = dt._ret;
    _pahitidx = dt._pahitidx;
    _sthitidx = dt._sthitidx;
  }

  TrkExtTraj & TrkExtTraj::operator= (const TrkExtTraj & dt) {
    _pt = dt._pt;
    _ret = dt._ret;
    _pahitidx = dt._pahitidx;
    _sthitidx = dt._sthitidx;
    return (*this);
  }

  TrkExtTrajPoint  & TrkExtTraj::operator[] (int i) {
    static TrkExtTrajPoint dummy;
    if (i <0 || i >= int(size())) return dummy;
    return _pt[i];
  }
      
  TrkExtTrajPoint    TrkExtTraj::operator[] (int i) const  {
    if (i <0 || i >= int(size())) return TrkExtTrajPoint();
    return _pt[i];
  }

  void TrkExtTraj::push_back (const TrkExtTrajPoint & trajPoint) {
    _pt.push_back(trajPoint);
  }

  void TrkExtTraj::clear () {
    _pt.clear();
    _pahitidx.clear();
    _sthitidx.clear();
  }



} // end namespace mu2e


