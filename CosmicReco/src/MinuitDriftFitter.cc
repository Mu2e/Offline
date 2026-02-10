// Author : S Middleton
// Date : August 2019
// Purpose: calls  minuit fitting to cosmic track seed. Input is CosmicTrackSeed, can then derive
// parameters from CosmicTrack stored there.

// ROOT:
#include "Offline/CosmicReco/inc/MinuitDriftFitter.hh"
#include "Offline/CosmicReco/inc/PDFFit.hh"
#include "Math/Math.h"
#include "Math/VectorUtil.h"
#include "Offline/Mu2eUtilities/inc/ParametricFit.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "TMath.h"
#include "Offline/TrackerGeom/inc/Tracker.hh"

// For Drift:
#include <TObjString.h>
#include <TROOT.h>
#include <TSystem.h>

// Minuit
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnPrint.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>

using namespace mu2e;

namespace MinuitDriftFitter {
FitResult DoFit(int const& _diag, CosmicTrackSeed& tseed, StrawResponse const& srep,
                const Tracker* tracker, double const& max_doca, unsigned int const& minChits,
                int const& MaxLogL, double const& _gaussTres, double const& maxTres) {

  std::vector<double> errors(5, 0);
  std::vector<double> seed(5, 0);
  std::vector<double> newseed(5, 0);
  std::vector<double> newerrors(5, 0);
  FitResult FitResult;
  CosmicTrack cosmictrack = tseed._track;

  // Seed Gaussian PDF using seed fit parameters stored in track info
  seed[0] = tseed._track.FitEquation.Pos.X();
  seed[1] = tseed._track.FitEquation.Dir.X();
  seed[2] = tseed._track.FitEquation.Pos.Y();
  seed[3] = tseed._track.FitEquation.Dir.Y();
  seed[4] = tseed._t0.t0();

  // Seed errors = covarience of parameters in seed fit
  errors[0] = tseed._track.FitParams.Covarience.sigA0;
  errors[1] = tseed._track.FitParams.Covarience.sigA1;
  errors[2] = tseed._track.FitParams.Covarience.sigB0;
  errors[3] = tseed._track.FitParams.Covarience.sigB1;
  errors[4] = tseed._t0.t0Err();

  // Constrain to mean = 0 for 4 parameters (T0 might not be so...)
  std::vector<double> constraint_means(5, 0);
  std::vector<double> constraints(5, 0);

  GaussianPDFFit fit(tseed._straw_chits, srep, cosmictrack, constraint_means, constraints,
                     _gaussTres, 1, tracker);

  // Initiate Minuit Fit:
  ROOT::Minuit2::MnStrategy mnStrategy(2);
  ROOT::Minuit2::MnUserParameters params(seed, errors);
  ROOT::Minuit2::MnMigrad migrad(fit, params, mnStrategy);

  // Set Limits as tracker dimensions:
  migrad.SetLimits((signed)0, -10000, 10000);
  migrad.SetLimits((signed)1, -5, 5);
  migrad.SetLimits((signed)2, -5000, 5000);
  migrad.SetLimits((signed)3, -10, 10);
  migrad.Fix((unsigned)4);
  int maxfcn = MaxLogL;
  double tolerance = 1000;

  // Define Minimization method as "MIGRAD" (see minuit documentation)
  ROOT::Minuit2::FunctionMinimum min = migrad(maxfcn, tolerance);
  if (_diag > 1) {
    ROOT::Minuit2::MnPrint::SetGlobalLevel(3);
    ROOT::Minuit2::operator<<(std::cout, min);
  }

  // Will be the results of the fit routine:
  ROOT::Minuit2::MnUserParameters results = min.UserParameters();
  double minval = min.Fval();

  // Define name for parameters
  FitResult.bestfit = results.Params();
  FitResult.bestfiterrors = results.Errors();

  // Store Minuit Covarience if exists:
  if (min.HasValidCovariance()) {
    FitResult.bestfitcov = min.UserCovariance().Data();
  }

  // Name Parameters:
  FitResult.names.push_back("a0");
  FitResult.names.push_back("a1");
  FitResult.names.push_back("b0");
  FitResult.names.push_back("b1");
  FitResult.names.push_back("t0");
  FitResult.NLL = minval;

  for (size_t i = 0; i < FitResult.names.size(); i++) {
    if ((!isnan(FitResult.NLL) or FitResult.NLL != 0 or FitResult.NLL < MaxLogL) and
        !isnan(FitResult.bestfit[i])) {
      cosmictrack.minuit_converged = true;
    }
  }

  if (_diag > 1) {
    for (size_t i = 0; i < FitResult.names.size(); i++) {
      std::cout << i << FitResult.names[i] << " : " << FitResult.bestfit[i] << " +- "
                << FitResult.bestfiterrors[i] << std::endl;
      if (FitResult.bestfitcov.size() != 0)
        std::cout << "cov " << FitResult.bestfitcov[i] << std::endl;
    }
  }

  // Cut on Gaussian results remove "bad" hits
  ComboHitCollection passed_hits;

  for (size_t i = 0; i < tseed._straw_chits.size(); i++) {
    double gauss_end_doca =
        fit.calculate_DOCA(tseed._straw_chits[i], FitResult.bestfit[0], FitResult.bestfit[1],
                           FitResult.bestfit[2], FitResult.bestfit[3], tracker);
    double gauss_end_time_residual = fit.TimeResidual(gauss_end_doca, srep, FitResult.bestfit[4],
                                                      tseed._straw_chits[i], tracker);
    FitResult.GaussianEndTimeResiduals.push_back(gauss_end_time_residual);
    FitResult.GaussianEndDOCAs.push_back(gauss_end_doca);
    if (gauss_end_doca < max_doca and gauss_end_time_residual < maxTres) {
      passed_hits.push_back(tseed._straw_chits[i]);
    } else {
      tseed._straw_chits[i]._flag.merge(StrawHitFlag::outlier);
    }
  }

  if (cosmictrack.minuit_converged and passed_hits.size() > minChits) {

    // Now run Full Fit:
    newseed[0] = FitResult.bestfit[0];
    newseed[1] = FitResult.bestfit[1];
    newseed[2] = FitResult.bestfit[2];
    newseed[3] = FitResult.bestfit[3];
    newseed[4] = FitResult.bestfit[4];

    newerrors[0] = FitResult.bestfiterrors[0];
    newerrors[1] = FitResult.bestfiterrors[1];
    newerrors[2] = FitResult.bestfiterrors[2];
    newerrors[3] = FitResult.bestfiterrors[3];
    newerrors[4] = tseed._t0.t0Err();
    FullDriftFit fulldriftfit(passed_hits, srep, cosmictrack, constraint_means, constraints,
                              _gaussTres, 1, tracker);

    ROOT::Minuit2::MnUserParameters newparams(newseed, newerrors);
    ROOT::Minuit2::MnMigrad newmigrad(fulldriftfit, newparams, mnStrategy);

    // Set Limits as tracker dimensions:
    newmigrad.SetLimits((signed)0, -10000, 10000);
    newmigrad.SetLimits((signed)1, -5, 5);
    newmigrad.SetLimits((signed)2, -5000, 5000);
    newmigrad.SetLimits((signed)3, -10, 10);
    newmigrad.Fix((unsigned)4);

    // Define Minimization method as "MIGRAD" (see minuit documentation)
    min = newmigrad(MaxLogL, tolerance);

    // Will be the results of the fit routine:
    results = min.UserParameters();
    minval = min.Fval();

    // Define name for parameters
    FitResult.bestfit = results.Params();
    FitResult.bestfiterrors = results.Errors();

    for (size_t i = 0; i < passed_hits.size(); i++) {

      double start_doca =
          fit.calculate_DOCA(passed_hits[i], seed[0], seed[1], seed[2], seed[3], tracker);
      double start_time_residual =
          fit.TimeResidual(start_doca, srep, seed[4], tseed._straw_chits[i], tracker);
      double end_doca =
          fit.calculate_DOCA(passed_hits[i], FitResult.bestfit[0], FitResult.bestfit[1],
                             FitResult.bestfit[2], FitResult.bestfit[3], tracker);
      double ambig = fit.calculate_ambig(passed_hits[i], FitResult.bestfit[0], FitResult.bestfit[1],
                                         FitResult.bestfit[2], FitResult.bestfit[3], tracker);
      double end_time_residual =
          fit.TimeResidual(end_doca, srep, FitResult.bestfit[4], passed_hits[i], tracker);

      FitResult.StartDOCAs.push_back(start_doca);
      FitResult.StartTimeResiduals.push_back(start_time_residual);
      FitResult.FullFitEndTimeResiduals.push_back(end_time_residual);
      FitResult.FullFitEndDOCAs.push_back(end_doca);
      FitResult.RecoAmbigs.push_back(ambig);
    }

    // delete array list to avoid memory leaks:
    fulldriftfit.DeleteArrays();
  }
  return FitResult;
}

void DoDriftTimeFit(
    std::vector<double> & pars,
    std::vector<double> & errors,
    std::vector<double> & cov_out,
    bool & minuit_converged,
    GaussianDriftFit &fit,
    double driftres,
    int diag, double mntolerance, double mnprecision) {

  // Initiate Minuit Fit:
  ROOT::Minuit2::MnStrategy mnStrategy(2);
  ROOT::Minuit2::MnUserParameters params(pars, errors);
  ROOT::Minuit2::MnMigrad migrad(fit, params, mnStrategy);

  if (mnprecision > 0) {
    migrad.SetPrecision(mnprecision);
  }

  // Do first fit stage with fixed drift res
  // and minimal t0
  fit.setFixedT0(true);
  fit.setFixedDriftRes(true,driftres);
  migrad.Fix(4);
  ROOT::Minuit2::FunctionMinimum temp_min = migrad(0, mntolerance);
  ROOT::Minuit2::MnUserParameters const& temp_results = temp_min.UserParameters();
  for (size_t i=0;i<4;i++){
    pars[i] = temp_results.Params()[i];
    errors[i] = temp_results.Errors()[i];
  }
  XYZVectorF ft0pos(pars[0], 0, pars[1]);
  XYZVectorF ft0dir(pars[2], -1, pars[3]);
  ft0dir = ft0dir.unit();
  pars[4] = fit.averageT0(pars);
  fit.setFixedT0(false);
  fit.setFixedDriftRes(false);
  migrad.Release(4);

  // Define Minimization method as "MIGRAD" (see minuit documentation)
  ROOT::Minuit2::FunctionMinimum min = migrad(0, mntolerance);
  if (diag > 1) {
    ROOT::Minuit2::MnPrint::SetGlobalLevel(3);
    ROOT::Minuit2::operator<<(std::cout, min);
  } else {
    ROOT::Minuit2::MnPrint::SetGlobalLevel(0);
  }

  // Will be the results of the fit routine:
  ROOT::Minuit2::MnUserParameters const& results = min.UserParameters();

  minuit_converged = min.IsValid();
  pars = results.Params();
  errors = results.Errors();

  if (min.HasValidCovariance()) {
    cov_out = min.UserCovariance().Data();
  } else {
    cov_out = std::vector<double>(15, 0);
  }
}

void DoDriftTimeFit(int const& diag, CosmicTrackSeed& tseed, StrawResponse const& srep,
                    const Tracker* tracker, double driftres, double mntolerance, double mnprecision) {

  auto dir = tseed._track.FitEquation.Dir;
  auto intercept = tseed._track.FitEquation.Pos;
  dir /= -1 * dir.y();
  intercept -= dir * intercept.y() / dir.y();

  // now gaussian fit, transverse distance only
  std::vector<double> errors(5, 0);
  std::vector<double> pars(5, 0);

  pars[0] = intercept.x();
  pars[1] = intercept.z();
  pars[2] = dir.x();
  pars[3] = dir.z();
  pars[4] = tseed._t0._t0;
  errors[0] = tseed._track.FitParams.Covarience.sigA0;
  errors[1] = tseed._track.FitParams.Covarience.sigB0;
  errors[2] = tseed._track.FitParams.Covarience.sigA1;
  errors[3] = tseed._track.FitParams.Covarience.sigB1;
  errors[4] = tseed._t0.t0Err();

  // Define the PDF used by Minuit:
  GaussianDriftFit fit(tseed._straw_chits, srep, tracker);
  DoDriftTimeFit(pars, errors, tseed._track.MinuitParams.cov,
    tseed._track.minuit_converged, fit, driftres,
    diag, mntolerance, mnprecision);

  tseed._track.MinuitParams.A0 = pars[0];
  tseed._track.MinuitParams.B0 = pars[1];
  tseed._track.MinuitParams.A1 = pars[2];
  tseed._track.MinuitParams.B1 = pars[3];
  tseed._track.MinuitParams.T0 = pars[4];
  tseed._track.MinuitParams.deltaA0 = errors[0];
  tseed._track.MinuitParams.deltaB0 = errors[1];
  tseed._track.MinuitParams.deltaA1 = errors[2];
  tseed._track.MinuitParams.deltaB1 = errors[3];
  tseed._track.MinuitParams.deltaT0 = errors[4];
  tseed._t0._t0 = tseed._track.MinuitParams.T0;
  tseed._t0._t0err = tseed._track.MinuitParams.deltaT0;

  XYZVectorF X(1, 0, 0);
  XYZVectorF Y(0, 1, 0);
  XYZVectorF Z(0, 0, 1);

  TrackAxes XYZ(X, Y, Z);
  tseed._track.MinuitCoordSystem = XYZ;
  tseed._track.MinuitEquation.Pos =
      XYZVectorF(tseed._track.MinuitParams.A0, 0, tseed._track.MinuitParams.B0);
  tseed._track.MinuitEquation.Dir =
      XYZVectorF(tseed._track.MinuitParams.A1, -1, tseed._track.MinuitParams.B1);

  for (size_t i = 0; i < tseed._straw_chits.size(); i++) {
    Straw const& straw = tracker->getStraw(tseed._straw_chits[i].strawId());
    TwoLinePCA pca(straw.getMidPoint(), straw.getDirection(),
                   GenVector::Hep3Vec(tseed._track.MinuitEquation.Pos),
                   GenVector::Hep3Vec(tseed._track.MinuitEquation.Dir));
    if (pca.dca() > 2.5) {
      tseed._straw_chits[i]._flag.merge(StrawHitFlag::outlier);
    }
  }
}

} // namespace MinuitDriftFitter
