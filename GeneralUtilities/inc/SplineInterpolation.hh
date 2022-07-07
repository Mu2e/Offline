#ifndef GeneralUtilities_SplineInterpolation_hh
#define GeneralUtilities_SplineInterpolation_hh

#include <vector>

// this spline will interpolate between xvals[0] and xvals[xvals.size()-1]
// making sure that it agrees with yvals at each point
// the spline arrays size will be xvals.size()-1
// calculates cubic hermite spline, optionally enforcing montonicity

namespace mu2e {

  class SplineInterpolation{

  public:
    SplineInterpolation(std::vector<double> const& xvals, std::vector<double> const& yvals,bool extrapolate, bool monotonic);
    SplineInterpolation(std::vector<double> const& xvals, std::vector<double> const& splineA,
        std::vector<double> const& splineB, std::vector<double> const& splineC,
        std::vector<double> const& splineD, bool extrapolate);
    ~SplineInterpolation(){};

    void getBin(double xval, int &ibin, double &t) const;

    double interpolate(int ibin, double t) const;
    double interpolate(double xval) const;
    double derivative(int ibin, double t) const;
    double derivative(double xval) const;

    std::vector<double> const& getSplineA(){return _splineA;}
    std::vector<double> const& getSplineB(){return _splineB;}
    std::vector<double> const& getSplineC(){return _splineC;}
    std::vector<double> const& getSplineD(){return _splineD;}

  private:
    std::vector<double> _xvals;
    double _deltax;
    std::vector<double> _splineA,_splineB,_splineC,_splineD;
    bool _extrapolate;
  };
}
#endif
