#ifndef _COSMIC_RECO_PDFFit_HH
#define _COSMIC_RECO_PDFFit_HH

#include "Offline/DataProducts/inc/GenVector.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh"
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"

// Tracker Details:
#include "Offline/TrackerConditions/inc/StrawDrift.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

// ROOT
#include "TF1.h"
#include "TH1F.h"
#include "TMath.h"

// Minuit
#include <Minuit2/FCNBase.h>

using namespace mu2e;

class GaussianPDFFit : public ROOT::Minuit2::FCNBase {
public:
  ComboHitCollection chits;
  StrawResponse const& srep;
  CosmicTrack track;
  std::vector<double> docas;

  std::vector<double> constraint_means;
  std::vector<double> constraints;
  double sigma_t;
  int k;

  const Tracker* tracker;
  std::vector<double> DOCAs;
  std::vector<double> TimeResiduals;

  int nparams = 5;

  GaussianPDFFit(ComboHitCollection const& _chits, StrawResponse const& _srep,
                 CosmicTrack const& _track, std::vector<double> const& _constraint_means,
                 std::vector<double> const& _constraints, double const& _sigma_t, int const& _k,
                 const Tracker* _tracker) :
      chits(_chits),
      srep(_srep), track(_track), constraint_means(_constraint_means), constraints(_constraints),
      sigma_t(_sigma_t), k(_k), tracker(_tracker){};

  double Up() const { return 0.5; };
  double operator()(const std::vector<double>& x) const;
  double TimeResidual(double const& doca, StrawResponse const& srep, double const& t0, ComboHit const& hit,
                      const Tracker* tracker) const;
  double calculate_DOCA(ComboHit const& chit, double const& a0, double const& a1, double const& b0,
                        double const& b1, const Tracker* tracker) const;
  double calculate_ambig(ComboHit const& chit, double const& a0, double const& a1, double const& b0,
                         double const& b1, const Tracker* tracker) const;
};

class FullDriftFit : public GaussianPDFFit {
public:
  FullDriftFit(ComboHitCollection const& _chits, StrawResponse const& _srep,
               CosmicTrack const& _track, std::vector<double> const& _constraint_means,
               std::vector<double> const& _constraints, double const& sigma_t, int const& _k,
               const Tracker* _tracker);
  int Factorial(int const& k);
  void CalculateFullPDF();
  double InterpolatePDF(double const& time_residual, double const& sigma, double const& tau) const;
  void DeleteArrays() const;
  double Min_t;
  double Min_tau, Min_s;
  double Max_t, Max_tau, Max_s;
  double delta_T, delta_Tau, delta_S;
  double *pdf_sigmas, *pdf_taus, *pdf_times;
  double* pdf;
  int k;
  double operator()(const std::vector<double>& x) const;
};

class GaussianDriftFit : public ROOT::Minuit2::FCNBase {
public:
  ComboHitCollection shs;
  StrawResponse const& srep;
  const Tracker* tracker;

  int excludeHit;

  bool fixedT0;
  bool fixedDriftRes;
  double driftRes;

  GaussianDriftFit(ComboHitCollection const& _shs, StrawResponse const& _srep,
                   const Tracker* _tracker) :
      shs(_shs),
      srep(_srep), tracker(_tracker), excludeHit(-1),
      fixedT0(false), fixedDriftRes(false), driftRes(0) {};
  // this tells Minuit to scale variances as if operator() returns a chi2 instead of a log
  // likelihood
  double Up() const { return 1.0; };
  double operator()(const std::vector<double>& x) const;

  void setExcludeHit(int const& hitIdx) {
    excludeHit = hitIdx;
  }
  void setFixedT0(bool fix) { fixedT0 = fix;}
  void setFixedDriftRes(bool fix, double dr=10) { fixedDriftRes = fix; driftRes = dr; }

  double averageT0(const std::vector<double> &x) const;

  double DOCAresidual(ComboHit const& sh, CosmicTrackSeed const& tseed) const {
    std::vector<double> x = {tseed._track.MinuitParams.A0, tseed._track.MinuitParams.B0,
                             tseed._track.MinuitParams.A1, tseed._track.MinuitParams.B1,
                             tseed._track.MinuitParams.T0};
    return DOCAresidual(sh, x);
  }
  double TimeResidual(ComboHit const& sh, CosmicTrackSeed const& tseed) const {
    std::vector<double> x = {tseed._track.MinuitParams.A0, tseed._track.MinuitParams.B0,
                             tseed._track.MinuitParams.A1, tseed._track.MinuitParams.B1,
                             tseed._track.MinuitParams.T0};
    return TimeResidual(sh, x);
  }
  double DOCAresidualError(ComboHit const& sh, CosmicTrackSeed const& tseed) const {
    std::vector<double> x = {tseed._track.MinuitParams.A0, tseed._track.MinuitParams.B0,
                             tseed._track.MinuitParams.A1, tseed._track.MinuitParams.B1,
                             tseed._track.MinuitParams.T0};
    return DOCAresidualError(sh, x, tseed._track.MinuitParams.cov);
  }

  int HitAmbiguity(ComboHit const& sh, CosmicTrackSeed const& tseed) const {
    std::vector<double> x {tseed._track.MinuitParams.A0, tseed._track.MinuitParams.B0,
                             tseed._track.MinuitParams.A1, tseed._track.MinuitParams.B1,
                             tseed._track.MinuitParams.T0};
    return HitAmbiguity(sh, x);
  }

  double reduced_chisq(const std::vector<double>& x);

  int HitAmbiguity(ComboHit const& sh, const std::vector<double>& x) const;

  double DOCAresidual(ComboHit const& sh, const std::vector<double>& x) const;

  double TimeResidual(ComboHit const& sh, const std::vector<double>& x) const;

  double DOCAresidualError(ComboHit const& sh, const std::vector<double>& x,
                           const std::vector<double>& cov) const;
};

#endif
