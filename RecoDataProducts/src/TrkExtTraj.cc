//
//  $Id: TrkExtTraj.cc,v 1.5 2013/02/07 02:09:47 mjlee Exp $
//  $Author: mjlee $
//  $Date: 2013/02/07 02:09:47 $
//
//  Original author MyeongJae Lee
//

#include "RecoDataProducts/inc/TrkExtTraj.hh"
#include "RecoDataProducts/inc/TrkExtTrajPoint.hh"
#include "CLHEP/Vector/ThreeVector.h"
#include "cetlib_except/exception.h"
#include <utility>
#include <cmath>

using namespace CLHEP;
using namespace std;

namespace mu2e {

  // Constructors
  TrkExtTraj::TrkExtTraj() {
    _pt.clear();
    _pahitidx.clear();
    _sthitidx.clear();
    _ptidx_pa.clear();
    _ptidx_st.clear();
    _deltap_pa.clear();
    _deltap_st.clear();
    _ret = -1;
    _hepid = 0;
  }

  TrkExtTraj::TrkExtTraj(const TrkExtTrajPoint & trajPoint) {
    _pt.clear();
    _pt.push_back(trajPoint);
    _pahitidx.clear();
    _sthitidx.clear();
    _ptidx_pa.clear();
    _ptidx_st.clear();
    _deltap_pa.clear();
    _deltap_st.clear();
    _ret = -1;
    _hepid = 0;
  }

  TrkExtTraj::TrkExtTraj(const TrkExtTraj & dt) {
    _pt = dt._pt;
    _ret = dt._ret;
    _hepid = dt._hepid;
    _pahitidx = dt._pahitidx;
    _sthitidx = dt._sthitidx;
    _ptidx_pa = dt._ptidx_pa;
    _ptidx_st = dt._ptidx_st;;
    _deltap_pa = dt._deltap_pa;
    _deltap_st = dt._deltap_st;
  }

  // operator overloading
  TrkExtTraj & TrkExtTraj::operator= (const TrkExtTraj & dt) {
    _pt = dt._pt;
    _ret = dt._ret;
    _hepid = dt._hepid;
    _pahitidx = dt._pahitidx;
    _sthitidx = dt._sthitidx;
    _ptidx_pa = dt._ptidx_pa;
    _ptidx_st = dt._ptidx_st;;
    _deltap_pa = dt._deltap_pa;
    _deltap_st = dt._deltap_st;
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
    _ptidx_pa.clear();
    _ptidx_st.clear();
    _deltap_pa.clear();
    _deltap_st.clear();
    _hepid = 0;
    _ret = -1;
  }

  // PA/ST related
  void TrkExtTraj::makePASTHitTable () {
    // see the note at header file for convention of first and second
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
    throw cet::exception("DATAPRODUCT")
      <<  "TrkExtTraj Error: Cannot find matching entry " << idx << ". PA/ST data now invalid. Report this problem to original author with reproducible example." << endl;
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


  // interpolation related
  unsigned int TrkExtTraj::findNeighborsZ (double z, std::vector<unsigned int> & idx, unsigned int idx1, unsigned int idx2) const{
    idx.clear();
    if (size() <=2) return 0;  // no interpolation when N(ext data) <=2
    bool p;
    unsigned int scanmin, scanmax;
    if (idx1 < idx2) {
      scanmin = idx1;
      scanmax = idx2;
    }
    else if (idx1 > idx2) {
      scanmin = idx2;
      scanmax = idx1;
    }
    else {
      scanmin = 1;
      scanmax = size()-1;
    }
    p = ( (_pt[scanmin].z() > z) ? true : false);
    for (unsigned int i = scanmin+1 ; i <=scanmax ; ++i) {
      bool n = ( (_pt[i].z() > z) ? true : false);
      if (n != p) {
        if ( fabs(_pt[i-1].z() - z) < fabs(_pt[i].z() - z) ) idx.push_back(i-1);
        else idx.push_back(i);
        p = n;
      }
    }
    if (idx.size()<=0) return 0;
    else return idx.size();
  }


  std::vector<TrkExtTrajPoint> TrkExtTraj::getPointsAtZ (double z, unsigned int idx1, unsigned int idx2) const {
    std::vector<TrkExtTrajPoint> points;
    points.clear();
    vector<unsigned int> idx;
    unsigned int ret = findNeighborsZ (z, idx, idx1, idx2);

    if (ret <=0) { ; }
/* two point interpolation
    else if (ret == 0) {
      unsigned int id = idx[0];
      double x = interpolate2(z, _pt[id].z(), _pt[id+1].z(), _pt[id].x(), _pt[id+1].x());
      double y = interpolate2(z, _pt[id].z(), _pt[id+1].z(), _pt[id].y(), _pt[id+1].y());
      double px = interpolate2(z, _pt[id].z(), _pt[id+1].z(), _pt[id].px(), _pt[id+1].px());
      double py = interpolate2(z, _pt[id].z(), _pt[id+1].z(), _pt[id].py(), _pt[id+1].py());
      double pz = interpolate2(z, _pt[id].z(), _pt[id+1].z(), _pt[id].pz(), _pt[id+1].pz());
      double fl = interpolate2(z, _pt[id].z(), _pt[id+1].z(), _pt[id].flightLength(), _pt[id+1].flightLength());
      double ft = interpolate2(z, _pt[id].z(), _pt[id+1].z(), _pt[id].flightTime(), _pt[id+1].flightTime());
      int volid;
      if (fabs(_pt[id].z() - z) < fabs(_pt[id+1].z() - z)) volid = _pt[id].volumeId();
      else volid = _pt[id+1].volumeId();
      Hep3Vector pos(x, y, z);
      Hep3Vector mom(px,py,pz);
      HepMatrix cov;
      TrkExtTrajPoint point(idx[0], pos, mom, cov, volid, fl, ft);
      points.push_back(point);
    }
*/
    else {
      for (unsigned int i = 0 ; i < idx.size() ; ++i) {
        unsigned int id = idx[i];
        double x = interpolate3(z, _pt[id-1].z(), _pt[id].z(), _pt[id+1].z(), _pt[id-1].x(), _pt[id].x(), _pt[id+1].x());
        double y = interpolate3(z, _pt[id-1].z(), _pt[id].z(), _pt[id+1].z(), _pt[id-1].y(), _pt[id].y(), _pt[id+1].y());
        double px = interpolate3(z, _pt[id-1].z(), _pt[id].z(), _pt[id+1].z(), _pt[id-1].px(), _pt[id].px(), _pt[id+1].px());
        double py = interpolate3(z, _pt[id-1].z(), _pt[id].z(), _pt[id+1].z(), _pt[id-1].py(), _pt[id].py(), _pt[id+1].py());
        double pz = interpolate3(z, _pt[id-1].z(), _pt[id].z(), _pt[id+1].z(), _pt[id-1].pz(), _pt[id].pz(), _pt[id+1].pz());
        double fl = interpolate3(z, _pt[id-1].z(), _pt[id].z(), _pt[id+1].z(), _pt[id-1].flightLength(), _pt[id].flightLength(), _pt[id+1].flightLength());
        double ft = interpolate3(z, _pt[id-1].z(), _pt[id].z(), _pt[id+1].z(), _pt[id-1].flightTime(), _pt[id].flightTime(), _pt[id+1].flightTime());
        int volid = _pt[id].volumeId();
        Hep3Vector pos(x, y, z);
        Hep3Vector mom(px,py,pz);
        HepMatrix cov(6,6,0);
        HepMatrix const & cov1 = _pt[id-1].covariance();
        HepMatrix const & cov2 = _pt[id].covariance();
        HepMatrix const & cov3 = _pt[id+1].covariance();
        for (unsigned int j = 1 ; j <=6 ; ++j) {
          for (unsigned int k = j ; k <=6 ; ++k) {
            cov(j,k) = cov(k,j) = interpolate3(z, _pt[id-1].z(), _pt[id].z(), _pt[id+1].z(), cov1(j,k), cov2(j,k), cov3(j,k));
          }
        }
        TrkExtTrajPoint point(idx[i], pos, mom, cov, volid, fl, ft);
        points.push_back(point);
      }
    }
    return points;
  }

/*  double TrkExtTraj::interpolate2 (double z, double x1, double x2, double y1, double y2) const {
    if (x2 != x1) {
      return (y2-y1)/(x2-x1)*(z-x1)+y1;
    }
    else {
      return 0;
    }
  }*/

  double TrkExtTraj::interpolate3 (double z, double x1, double x2, double x3, double y1, double y2, double y3) const {
    double x1_x2 = x1 - x2;
    double x2_x3 = x2 - x3;
    double x3_x1 = x3 - x1;
    if (x1_x2 !=0 && x2_x3 != 0 && x3_x1 != 0) {
      return y1*(z-x2)*(x3-z)/x1_x2/x3_x1 + y2*(x1-z)*(z-x3)/x1_x2/x2_x3 + y3*(z-x1)*(x2-z)/x3_x1/x2_x3;
    }
    else {
      return 0;
    }
  }


  const TrkExtTrajPoint & TrkExtTraj::getFirstPAHit (unsigned int i) const {
    // always returns final hit inside PA in time
    int first = _ptidx_pa[i].first;
    int second = _ptidx_pa[i].second;

    int nsteps = std::abs(second - first);
    if (nsteps ==0) {
      return _pt[first];
    }
    else {
      if (_pt[first].flightTime() > _pt[second].flightTime())
        return _pt[first];
      else
        return _pt[second];
    }
  }

  const TrkExtTrajPoint & TrkExtTraj::getFirstSTHit (unsigned int i) const {
    // always returns final hit inside ST in time
    int first = _ptidx_st[i].first;
    int second = _ptidx_st[i].second;

    int nsteps = std::abs(second - first);
    if (nsteps ==0) {
      return _pt[first];
    }
    else {
      if (_pt[first].flightTime() > _pt[second].flightTime())
        return _pt[first];
      else
        return _pt[second];
    }
  }


  TrkExtTrajPoint TrkExtTraj::getMeanPAHit (unsigned int i) const {
    // returns interpolated center inside PA in z
    unsigned int first = _ptidx_pa[i].first;
    unsigned int second = _ptidx_pa[i].second;

    if (first == second) return _pt[first];

    double zc = (_pt[first].z() + _pt[second].z()) * 0.5;

    vector<TrkExtTrajPoint> points = getPointsAtZ(zc, first, second);
    if (points.size() <=0) {
      throw cet::exception("DATAPRODUCT")
        << "TrkExtTraj Error : cannot find getMeanPaHit. Error in interpolation! Report to the author." << endl;
      TrkExtTrajPoint point;
      return point;
    }
    return points[0];
  }

  TrkExtTrajPoint TrkExtTraj::getMeanSTHit (unsigned int i) const {
    // returns interpolated center inside ST in z
    unsigned int first = _ptidx_st[i].first;
    unsigned int second = _ptidx_st[i].second;

    if (first == second) return _pt[first];

    double zc = (_pt[first].z() + _pt[second].z()) * 0.5;

    vector<TrkExtTrajPoint> points = getPointsAtZ(zc, first, second);
    if (points.size() <=0) {
      throw cet::exception("DATAPRODUCT")
        << "TrkExtTraj Error : cannot find getMeanPaHit. Error in interpolation! Report to the author." << endl;
      TrkExtTrajPoint point;
      return point;
    }
    return points[0];
  }



} // end namespace mu2e
