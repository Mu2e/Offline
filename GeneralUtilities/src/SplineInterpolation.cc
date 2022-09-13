#include "Offline/GeneralUtilities/inc/SplineInterpolation.hh"

#include <cmath>
#include <algorithm>

using namespace std;

namespace mu2e {
  SplineInterpolation::SplineInterpolation(std::vector<double> const& xvals, std::vector<double> const& yvals, bool extrapolate, bool monotonic) :
    _xvals(xvals), _extrapolate(extrapolate)
  {
    _deltax = (xvals.back()-xvals.front())/(xvals.size()-1);
    _splineA.resize(xvals.size()-1);
    _splineB.resize(xvals.size()-1);
    _splineC.resize(xvals.size()-1);
    _splineD.resize(xvals.size()-1);

    for (size_t i=0;i<xvals.size()-1;i++){
      double p0,p1,p2,p3,m1,m2;
      if (i == 0){
        p0 = (yvals[i] - (yvals[i+1]-yvals[i]));
      }else{
        p0 = yvals[i-1];
      }
      if (i == xvals.size()-2){
        p3 = (yvals[i+1] + (yvals[i+1]-yvals[i]));
      }else{
        p3 = yvals[i+2];
      }
      p1 = yvals[i];
      p2 = yvals[i+1];
      m1 = (p2-p0)/2.;
      m2 = (p3-p1)/2.;

      if (monotonic){
        if ((p2-p1)*(p1-p0) <= 0){
          m1 = 0;
        }
        if ((p3-p2)*(p2-p1) <= 0){
          m2 = 0;
        }
        if (m1 > 0){
          m1 = min(max(0.0,m1),3*min(p2-p1,p1-p0));
        }else if (m1 < 0){
          m1 = max(min(0.0,m1),3*max(p2-p1,p1-p0));
        }
        if (m2 > 0){
          m2 = min(max(0.0,m2),3*min(p3-p2,p2-p1));
        }else if (m2 < 0){
          m2 = max(min(0.0,m2),3*max(p3-p2,p2-p1));
        }
      }

      _splineD[i] = p1;
      _splineC[i] = m1;
      _splineB[i] = -3*p1 + 3*p2 - 2*m1 - m2;
      _splineA[i] = 2*p1 + m1 - 2*p2 + m2;
    }
  }

  SplineInterpolation::SplineInterpolation(std::vector<double> const& xvals, std::vector<double> const& splineA,
        std::vector<double> const& splineB, std::vector<double> const& splineC,
        std::vector<double> const& splineD, bool extrapolate) :
    _xvals(xvals), _splineA(splineA), _splineB(splineB), _splineC(splineC), _splineD(splineD), _extrapolate(extrapolate)
  {
    if (xvals.size() > 1)
      _deltax = (xvals.back()-xvals.front())/(xvals.size()-1);
  }

  void SplineInterpolation::getBin(double xval, int &ibin, double &t) const {
    int imax = int(_xvals.size()-2);
    ibin = min(imax,max(0,int(floor((xval-_xvals.front())/_deltax))));
    t = (xval - _xvals[ibin])/_deltax;
  }

  double SplineInterpolation::interpolate(int ibin, double t) const {
    if (t >= 1){
      if (_extrapolate){
        return ((3*(t-1) + 1)*_splineA[ibin] + (2*(t-1) + 1)*_splineB[ibin] + ((t-1) + 1)*_splineC[ibin]) + _splineD[ibin];
      }else{
        return _splineA[ibin] + _splineB[ibin] + _splineC[ibin] + _splineD[ibin];
      }
    }
    return _splineA[ibin]*pow(t,3) + _splineB[ibin]*pow(t,2) + _splineC[ibin]*t + _splineD[ibin];
  }

  double SplineInterpolation::interpolate(double xval) const {
    int ibin;
    double t;
    getBin(xval,ibin,t);
    return interpolate(ibin,t);
  }

  double SplineInterpolation::derivative(int ibin, double t) const {
    if (t >= 1){
      if (_extrapolate){
        return (3*_splineA[ibin]+2*_splineB[ibin]+_splineC[ibin])/_deltax;
      }else{
        return 0;
      }
    }
    return (3*_splineA[ibin]*pow(t,2) + 2*_splineB[ibin]*t + _splineC[ibin])/_deltax;
  }

  double SplineInterpolation::derivative(double xval) const {
    int ibin;
    double t;
    getBin(xval,ibin,t);
    return derivative(ibin,t);
  }

  /*
  double SplineInterpolation::interpolate(double xval) const {
    int imax = int(_xvals.size()-2);
    if (xval >= _xvals[_xvals.size()-1]){
      double t = (xval - _xvals[_xvals.size()-1])/_deltax;
      return ((3*t + 1)*_splineA[imax] + (2*t + 1)*_splineB[imax] + (t + 1)*_splineC[imax]) + _splineD[imax];
    }
    int ibin = min(imax,max(0,int(floor((xval-_xvals.front())/_deltax))));
    double t = (xval - (double) (ibin * _deltax))/_deltax;
    return _splineA[ibin]*pow(t,3) + _splineB[ibin]*pow(t,2) + _splineC[ibin]*t + _splineD[ibin];
  }

  double SplineInterpolation::derivative(double xval) const {
    int imax = int(_xvals.size()-2);
    if (xval >= _xvals[_xvals.size()-1]){
      return (3*_splineA[imax]+2*_splineB[imax]+_splineC[imax])/_deltax;
    }
    int ibin = min(imax,max(0,int(floor((xval-_xvals.front())/_deltax))));
    double t = (xval - (double) (ibin * _deltax))/_deltax;
    return (3*_splineA[ibin]*pow(t,2) + 2*_splineB[ibin]*t + _splineC[ibin])/_deltax;
  }
  */

}
