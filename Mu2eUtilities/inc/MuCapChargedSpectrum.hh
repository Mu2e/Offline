#ifndef Mu2eUtilities_MuCapChargedSpectrum_hh
#define Mu2eUtilities_MuCapChargedSpectrum_hh
//
// Kinetic energy spectrum of a charged particle (p, d, ...) emitted after muon capture, in
// absolute units: getWeight(T) is the yield per capture per MeV, so the integral is the rate.
// Usable with BinnedSpectrum (spectrumShape : muCapCharged).
//
// The form follows the TWIST parameterization (PRC 101, 035502 (2020), arXiv:1908.06902):
//
//   dN/dT = A exp(-gamma T) (1 + phi(p))     for p > p_min
//
// with phi a cubic B-spline on (p_min, p_u] and 0 above, p the momentum. Below p_min (T < T_min, no
// data) one of two onsets is used:
//
//   onset : { power : a }             k T^a exp(-b T), with k and b set so that the value AND the
//                                     slope are continuous at T_min; goes to zero smoothly at T -> 0
//   onset : { Tc, alpha, beta }       k (1 - Tc/T)^alpha exp(-beta T) for T > Tc, 0 below; k set by
//                                     continuity of the value only (first version)
//
// Optionally a hard component is added,
//
//   H(T) = h exp(-(T - Ton)/slope) / (1 + exp(-(T - Ton)/width)),
//
// with h chosen so that its integral from Ton to the endpoint is the configured yield (before the
// taper). The spectrum is zero above the kinematic endpoint. With taper : { width : w } the whole
// spectrum is multiplied by (1 - x^2)^2, x = (T - (endpoint - w))/w, over the last w MeV, so that it
// goes to zero smoothly, like (endpoint - T)^2, at the endpoint.
//
// Original author: Sophie Middleton, 2026
//
#include <vector>

namespace fhicl { class ParameterSet; }

namespace mu2e {

  class MuCapChargedSpectrum {
  public:
    struct Config {
      double mass = 0.;            // MeV/c^2
      double endpoint = 0.;        // maximum kinetic energy, MeV
      double A = 0.;               // 1/MeV
      double gamma = 0.;           // 1/MeV
      std::vector<double> knots;   // MeV/c, cubic: weights.size() + 4 entries
      std::vector<double> weights;
      double onsetPower = 0.;      // > 0: smooth C1 onset k T^power exp(-b T); else the (Tc, alpha, beta) onset
      double onsetTc = 0.;         // MeV
      double onsetAlpha = 0.;
      double onsetBeta = 0.;       // 1/MeV
      double taperWidth = 0.;      // MeV; 0: hard cut at the endpoint
      double hardYield = 0.;       // per capture, between hardTon and the endpoint
      double hardSlope = 1.;       // MeV
      double hardTon = 0.;         // MeV
      double hardWidth = 1.;       // MeV
    };

    explicit MuCapChargedSpectrum(const Config& config);
    // pset: the BinnedSpectrum parameter set (pdgId, endpoint, tail, spline, onset, [hard], [taper])
    explicit MuCapChargedSpectrum(const fhicl::ParameterSet& pset);

    double getWeight(double T) const; // T in MeV

    double endpoint() const { return _c.endpoint; }

    // default endpoint from the pset, for BinnedSpectrum before construction
    static double endpoint(const fhicl::ParameterSet& pset);

  private:
    static Config parse(const fhicl::ParameterSet& pset);
    void init();
    double measured(double p) const; // TWIST part, p >= p_min
    double phi(double p) const;
    double onset(double T) const;
    double hard(double T) const;
    double taper(double T) const;

    Config _c;
    double _Tmin = 0.;  // kinetic energy at the first knot
    double _onsetK = 0.;
    double _onsetB = 0.;  // smooth onset exponent, 1/MeV
    double _hardNorm = 0.;
  };

}

#endif /* Mu2eUtilities_MuCapChargedSpectrum_hh */
