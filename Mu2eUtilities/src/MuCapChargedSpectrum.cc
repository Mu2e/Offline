// See the header for the parameterization.
#include "Offline/Mu2eUtilities/inc/MuCapChargedSpectrum.hh"

#include "Offline/DataProducts/inc/PDGCode.hh"
#include "Offline/GlobalConstantsService/inc/GlobalConstantsHandle.hh"
#include "Offline/GlobalConstantsService/inc/ParticleDataList.hh"

#include "cetlib_except/exception.h"
#include "fhiclcpp/ParameterSet.h"

#include <cmath>

namespace {
  // Cox-de Boor recursion for basis function i of degree k, closed at the last knot
  double bspline(const std::vector<double>& t, size_t i, int k, double x) {
    if(k == 0) {
      if(t[i] <= x && x < t[i+1]) return 1.;
      if(x == t.back() && t[i+1] == t.back() && t[i] < t[i+1]) return 1.;
      return 0.;
    }
    double res = 0.;
    if(t[i+k] > t[i])     res += (x - t[i])/(t[i+k] - t[i]) * bspline(t, i, k-1, x);
    if(t[i+k+1] > t[i+1]) res += (t[i+k+1] - x)/(t[i+k+1] - t[i+1]) * bspline(t, i+1, k-1, x);
    return res;
  }
}

namespace mu2e {

  MuCapChargedSpectrum::MuCapChargedSpectrum(const Config& config) : _c(config) { init(); }

  MuCapChargedSpectrum::MuCapChargedSpectrum(const fhicl::ParameterSet& pset) : _c(parse(pset)) { init(); }

  double MuCapChargedSpectrum::endpoint(const fhicl::ParameterSet& pset) {
    return pset.get<double>("endpoint");
  }

  MuCapChargedSpectrum::Config MuCapChargedSpectrum::parse(const fhicl::ParameterSet& pset) {
    Config c;
    const auto pdgId = static_cast<PDGCode::type>(pset.get<int>("pdgId"));
    c.mass = GlobalConstantsHandle<ParticleDataList>()->particle(pdgId).mass();
    c.endpoint = pset.get<double>("endpoint");
    const auto tail = pset.get<fhicl::ParameterSet>("tail");
    c.A     = tail.get<double>("A");
    c.gamma = tail.get<double>("gamma");
    const auto spline = pset.get<fhicl::ParameterSet>("spline");
    c.knots   = spline.get<std::vector<double>>("knots");
    c.weights = spline.get<std::vector<double>>("weights");
    const auto onset = pset.get<fhicl::ParameterSet>("onset");
    c.onsetTc    = onset.get<double>("Tc");
    c.onsetAlpha = onset.get<double>("alpha");
    c.onsetBeta  = onset.get<double>("beta");
    const auto hard = pset.get<fhicl::ParameterSet>("hard", fhicl::ParameterSet());
    c.hardYield = hard.get<double>("yield", 0.);
    c.hardSlope = hard.get<double>("slope", 1.);
    c.hardTon   = hard.get<double>("Ton", 0.);
    c.hardWidth = hard.get<double>("width", 1.);
    return c;
  }

  void MuCapChargedSpectrum::init() {
    if(_c.knots.size() != _c.weights.size() + 4)
      throw cet::exception("BADCONFIG") << "MuCapChargedSpectrum: need weights.size()+4 spline knots, got "
                                        << _c.knots.size() << " for " << _c.weights.size() << " weights";
    if(!(_c.mass > 0.) || !(_c.endpoint > 0.) || !(_c.hardSlope > 0.) || !(_c.hardWidth > 0.))
      throw cet::exception("BADCONFIG") << "MuCapChargedSpectrum: mass, endpoint, hard slope and width must be > 0";

    const double pmin = _c.knots.front();
    _Tmin = std::sqrt(pmin*pmin + _c.mass*_c.mass) - _c.mass;
    if(_Tmin <= _c.onsetTc)
      throw cet::exception("BADCONFIG") << "MuCapChargedSpectrum: onset Tc " << _c.onsetTc
                                        << " must be below the first knot (T = " << _Tmin << " MeV)";
    _onsetK = 1.;
    _onsetK = measured(pmin)/onset(_Tmin); // continuity at p_min

    // normalize the hard component to its yield between Ton and the endpoint
    _hardNorm = 0.;
    if(_c.hardYield > 0. && _c.hardTon < _c.endpoint) {
      _hardNorm = 1.;
      const int n = 10000;
      const double dT = (_c.endpoint - _c.hardTon)/n;
      double sum = 0.;
      for(int i = 0; i < n; ++i) sum += hard(_c.hardTon + (i + 0.5)*dT);
      _hardNorm = _c.hardYield/(sum*dT);
    }
  }

  double MuCapChargedSpectrum::phi(double p) const {
    if(p < _c.knots.front() || p > _c.knots.back()) return 0.;
    double res = 0.;
    for(size_t l = 0; l < _c.weights.size(); ++l) res += _c.weights[l]*bspline(_c.knots, l, 3, p);
    return res;
  }

  double MuCapChargedSpectrum::measured(double p) const {
    const double T = std::sqrt(p*p + _c.mass*_c.mass) - _c.mass;
    return _c.A*std::exp(-_c.gamma*T)*(1. + phi(p));
  }

  double MuCapChargedSpectrum::onset(double T) const {
    if(T <= _c.onsetTc) return 0.;
    return _onsetK*std::pow(1. - _c.onsetTc/T, _c.onsetAlpha)*std::exp(-_c.onsetBeta*T);
  }

  double MuCapChargedSpectrum::hard(double T) const {
    if(_hardNorm == 0.) return 0.;
    const double x = (T - _c.hardTon);
    return _hardNorm*std::exp(-x/_c.hardSlope)/(1. + std::exp(-x/_c.hardWidth));
  }

  double MuCapChargedSpectrum::getWeight(double T) const {
    if(T <= 0. || T > _c.endpoint) return 0.;
    const double p = std::sqrt(T*(T + 2.*_c.mass));
    const double soft = (p < _c.knots.front()) ? onset(T) : measured(p);
    return soft + hard(T);
  }

}
