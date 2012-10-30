//
//  $Id: TrkExtTraj.cc,v 1.3 2012/10/30 22:05:27 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2012/10/30 22:05:27 $
//
//  Original author MyeongJae Lee
//
//


#include "RecoDataProducts/inc/TrkExtTraj.hh"
#include "RecoDataProducts/inc/TrkExtTrajPoint.hh"
#include <utility>

using namespace std;

namespace mu2e {

  // Constructors
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

  // operator overloading
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

  //vector-like accessor
  void TrkExtTraj::push_back (const TrkExtTrajPoint & trajPoint) {
    _pt.push_back(trajPoint);
  }

  void TrkExtTraj::clear () {
    _pt.clear();
    _pahitidx.clear();
    _sthitidx.clear();
  }

  // PA/ST related
  void TrkExtTraj::makePASTHitTable () {
    _ptidx_pa.clear();
    _ptidx_st.clear();
    _deltap_pa.clear();
    _deltap_st.clear();
    unsigned int ret = 0, first, second;
    double esum = 0;
    if (_pahitidx.size() >0) {
      for (unsigned int i = 0 ; i < _pahitidx.size() ; ++i) {
        ret = findPASTHit(_pahitidx[i].first, ret);
        first = ret;
        ret = findPASTHit(_pahitidx[i].second, ret+1);
        second = ret;
        _ptidx_pa.push_back(make_pair(first, second));
        esum = 0;
        if (second > first) esum = _pt[second].momentum().mag() - _pt[first].momentum().mag();
        _deltap_pa.push_back(esum);
      }
    }
    if (_sthitidx.size() >0) {
      for (unsigned int i = 0 ; i < _sthitidx.size() ; ++i) {
        ret = findPASTHit(_sthitidx[i].first, ret);
        first = ret;
        ret = findPASTHit(_sthitidx[i].second, ret+1);
        second = ret;
        _ptidx_st.push_back(make_pair(first, second));
        esum = 0;
        if (second > first) esum = _pt[second].momentum().mag() - _pt[first].momentum().mag();
        _deltap_st.push_back(esum);
      }
    }
  }

  unsigned int TrkExtTraj::findPASTHit (unsigned int idx, unsigned int start) {
    int iidx = (int)idx;
    for (unsigned int i = start ; i < _pt.size() ; ++i) {
      if (_pt[i].trajPointId() == iidx) return i;
    }
    for (unsigned int i = 0 ; i < start ; ++i) {
      if (_pt[i].trajPointId() == iidx) return i;
    }
    cout <<  "TrkExtTraj Error: Cannor find matching entry " << idx << " PA/ST data now invalid." << endl;
    return 0;
  }

  double TrkExtTraj::getDeltapPA () const{
    double esum = 0;
    for (unsigned int i  = 0 ; i < _deltap_pa.size() ; ++i) {
      esum += _deltap_pa[i];
    }
    return esum;
  }

  double TrkExtTraj::getDeltapST () const {
    double esum = 0;
    for (unsigned int i  = 0 ; i < _deltap_st.size() ; ++i) {
      esum += _deltap_st[i];
    }
    return esum;
  }





} // end namespace mu2e


