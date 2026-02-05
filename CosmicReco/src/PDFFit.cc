// Author: S Middleton
// Date: August 2019
// Purpose: PDF Functions for minuit fitting

// ROOT:
#include "Math/Math.h"
#include "Math/VectorUtil.h"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"
#include "TMath.h"
#include <TObjString.h>
#include <TROOT.h>
#include <TSystem.h>

// Tracker Details:
#include "Offline/TrackerConditions/inc/StrawDrift.hh"
#include "Offline/TrackerConditions/inc/StrawResponse.hh"
#include "Offline/TrackerGeom/inc/Tracker.hh"

// Utilities:
#include "Offline/CosmicReco/inc/DriftFitUtils.hh"
#include "Offline/CosmicReco/inc/PDFFit.hh"
#include "Offline/Mu2eUtilities/inc/TwoLinePCA.hh"

// Minuit
#include <Minuit2/FCNBase.h>
#include <Minuit2/FunctionMinimum.h>
#include <Minuit2/MnMigrad.h>
#include <Minuit2/MnMinos.h>
#include <Minuit2/MnStrategy.h>
#include <Minuit2/MnUserParameters.h>

using namespace mu2e;

// For Gaussian:
// const double drift_v = 0.065; //mm/ns TODO fix this in the code
const double sigma =
    1.0; // ns (for full fit only) , if gaussian the sigma is set to 24ns in the FCL parametets
const double tau = 10.7; // ns
// For Full fit:
const int N_tbins = 500;
const int N_taubins = 50;
const int N_sbins = 50;

float wireradius = 12.5 / 1000.; // 12.5 um in mm
float strawradius = 2.5;         // 2.5 mm in mm

FullDriftFit::FullDriftFit(ComboHitCollection const& _chits, StrawResponse const& _srep,
                           CosmicTrack const& _track, std::vector<double> const& _constraint_means,
                           std::vector<double> const& _constraints, double const& _sigma_t,
                           int const& _k, const Tracker* _tracker) :
    GaussianPDFFit(_chits, _srep, _track, _constraint_means, _constraints, _sigma_t, _k, _tracker) {
  pdf_times = new double[N_tbins];
  pdf_taus = new double[N_taubins];
  pdf_sigmas = new double[N_sbins];
  pdf = new double[N_tbins * N_taubins * N_sbins];

  Min_t = 0.1;
  Min_tau = 5;
  Min_s = -1;
  Max_t = 1000;
  Max_tau = 15.;
  Max_s = 30.;
  k = 1;

  delta_T = (Max_t - Min_t) / ((double)N_tbins);
  delta_S = (Max_s - Min_s) / ((double)N_sbins);
  delta_Tau = (Max_tau - Min_tau) / ((double)N_taubins);

  for (int i = 0; i < N_tbins; i++) {
    pdf_times[i] = Min_t + delta_T * i;
  }
  for (int i = 0; i < N_sbins; i++) {
    pdf_sigmas[i] = Min_s + delta_S * i;
  }
  for (int i = 0; i < N_taubins; i++) {
    pdf_taus[i] = Min_tau + delta_Tau * i;
  }
  this->CalculateFullPDF();
}

int FullDriftFit::Factorial(int const& k) {
  if (k == 0)
    return 1;
  int response = 1;
  for (int i = 1; i < k; i++) {
    response *= i;
  }
  return response;
}

void FullDriftFit::CalculateFullPDF() {

  for (int is = 0; is < N_sbins; is++) {
    double sigma = this->pdf_sigmas[is];
    for (int it0 = 0; it0 < N_tbins; it0++) {
      double time_gaus = this->pdf_times[it0];
      double time_gaussian = 1.0 / sqrt(2 * TMath::Pi() * sigma * sigma) *
                             exp(-(time_gaus * time_gaus) / (2 * sigma * sigma));

      for (int itau = 0; itau < N_taubins; itau++) {
        double tau = this->pdf_taus[itau];
        for (int it1 = 0; it1 < N_tbins - it0; it1++) {
          double time_tau = this->delta_T * it1;
          double val_tau = pow(1 / tau, k) * pow(time_tau, k - 1) * exp(-time_tau / tau) /
                           (double)Factorial(k - 1);
          this->pdf[is * N_taubins * N_tbins + itau * N_tbins + (it0 + it1)] +=
              time_gaussian * val_tau;
        }
      }
    }
  }

  for (int is = 0; is < N_sbins; is++) {
    for (int itau = 0; itau < N_taubins; itau++) {
      double total = 0;
      for (int it = 0; it < N_tbins; it++) {
        total += this->pdf[is * N_taubins * N_tbins + itau * N_tbins + it];
      }

      for (int it = 0; it < N_tbins; it++) {
        this->pdf[(is * N_taubins * N_tbins) + (itau * N_tbins) + it] /= total;
      }
    }
  }
}

void FullDriftFit::DeleteArrays() const {

  delete[] pdf_sigmas;
  delete[] pdf_taus;
  delete[] pdf_times;
  delete[] pdf;
}

double FullDriftFit::InterpolatePDF(double const& time_residual, double const& sigma,
                                    double const& tau) const {

  int bin_s = (sigma - this->Min_s) / (this->delta_S);
  if (bin_s >= N_sbins - 1) {
    bin_s = N_sbins - 2;
  }
  if (bin_s < 0) {
    bin_s = 0;
  }
  double s_d = (sigma - (bin_s * this->delta_S + this->Min_s)) / (this->delta_S);
  int bin_tau = (tau - this->Min_tau) / (this->delta_Tau);

  if (bin_tau >= N_taubins - 1) {
    bin_tau = N_taubins - 2;
  }

  double tau_d = (tau - (bin_tau * this->delta_Tau + this->Min_tau)) / (this->delta_Tau);
  int bin_t = (time_residual - this->Min_t) / (this->delta_T);

  if (bin_t >= N_tbins - 1) {
    bin_t = N_tbins - 2;
  }
  if (bin_t < 0) {
    bin_t = 0;
  }
  double t_d = (time_residual - (bin_t * this->delta_T + this->Min_t)) / (this->delta_T);

  double pdf_val000 =
      this->pdf[(bin_s + 0) * N_taubins * N_tbins + (bin_tau + 0) * N_tbins + (bin_t + 0)];

  double pdf_val100 =
      this->pdf[(bin_s + 1) * N_taubins * N_tbins + (bin_tau + 0) * N_tbins + (bin_t + 0)];

  double pdf_val001 =
      this->pdf[(bin_s + 0) * N_taubins * N_tbins + (bin_tau + 0) * N_tbins + (bin_t + 1)];

  double pdf_val101 =
      this->pdf[(bin_s + 1) * N_taubins * N_tbins + (bin_tau + 0) * N_tbins + (bin_t + 1)];

  double pdf_val010 =
      this->pdf[(bin_s + 0) * N_taubins * N_tbins + (bin_tau + 1) * N_tbins + (bin_t + 0)];

  double pdf_val110 =
      this->pdf[(bin_s + 1) * N_taubins * N_tbins + (bin_tau + 1) * N_tbins + (bin_t + 0)];

  double pdf_val011 =
      this->pdf[(bin_s + 0) * N_taubins * N_tbins + (bin_tau + 1) * N_tbins + (bin_t + 1)];

  double pdf_val111 =
      this->pdf[(bin_s + 1) * N_taubins * N_tbins + (bin_tau + 1) * N_tbins + (bin_t + 1)];

  double pdf_val00 = pdf_val000 * (1 - s_d) + pdf_val100 * s_d;
  double pdf_val01 = pdf_val001 * (1 - s_d) + pdf_val101 * s_d;
  double pdf_val10 = pdf_val010 * (1 - s_d) + pdf_val110 * s_d;
  double pdf_val11 = pdf_val011 * (1 - s_d) + pdf_val111 * s_d;
  double pdf_val0 = pdf_val00 * (1 - tau_d) + pdf_val10 * tau_d;
  double pdf_val1 = pdf_val01 * (1 - tau_d) + pdf_val11 * tau_d;
  double pdf_val = pdf_val0 * (1 - t_d) + pdf_val1 * t_d;

  return pdf_val;
}

// This 3 functions talk to the drift util:
double GaussianPDFFit::calculate_DOCA(ComboHit const& chit, double const& a0, double const& a1,
                                      double const& b0, double const& b1,
                                      const Tracker* tracker) const {
  double doca = DriftFitUtils::GetTestDOCA(chit, a0, a1, b0, b1, tracker);
  return (doca);
}

double GaussianPDFFit::calculate_ambig(ComboHit const& chit, double const& a0, double const& a1,
                                       double const& b0, double const& b1,
                                       const Tracker* tracker) const {
  double ambig = DriftFitUtils::GetAmbig(chit, a0, a1, b0, b1, tracker);
  return (ambig);
}

double GaussianPDFFit::TimeResidual(double const& doca, StrawResponse const& srep, double const& t0,
                                    ComboHit const& hit, const Tracker* tracker) const {
  double tres = DriftFitUtils::TimeResidual(doca, srep, t0, hit, tracker);
  return (tres);
}

double GaussianPDFFit::operator()(const std::vector<double>& x) const {

  double const& a0 = x[0];
  double const& a1 = x[1];
  double const& b0 = x[2];
  double const& b1 = x[3];
  double t0 = x[4];
  long double llike = 0;

  for (size_t i = 0; i < this->chits.size(); i++) {
    double doca = calculate_DOCA(this->chits[i], a0, a1, b0, b1, this->tracker);
    double time_residual = this->TimeResidual(doca, this->srep, t0, this->chits[i], this->tracker);
    double pdf_t = 1 / sqrt(2 * TMath::Pi() * this->sigma_t * this->sigma_t) *
                   exp(-(time_residual * time_residual) / (2 * this->sigma_t * this->sigma_t));

    llike -= log(pdf_t);
    t0 += time_residual / this->chits.size();
  }

  for (int i = 0; i < this->nparams; i++) {
    if (this->constraints[i] > 0) {
      llike += pow((x[i] - this->constraint_means[i]) / this->constraints[i], 2);
    }
  }

  return (double)llike;
}

double FullDriftFit::operator()(const std::vector<double>& x) const {

  if (sigma > Max_s) {
    return 1e10;
  }

  double const& a0 = x[0];
  double const& a1 = x[1];
  double const& b0 = x[2];
  double const& b1 = x[3];
  double t0 = x[4];
  long double llike = 0;

  for (size_t i = 0; i < this->chits.size(); i++) {
    double doca = calculate_DOCA(this->chits[i], a0, a1, b0, b1, this->tracker);
    double doca_penalty = 1;
    double time_residual = this->TimeResidual(doca, this->srep, t0, this->chits[i], this->tracker);

    if (time_residual > 40) {
      time_residual = 40;
    }

    double hypotenuse = sqrt(pow(doca, 2) + pow((tau * 0.0625), 2));

    double tau_eff = (hypotenuse / 0.0625) - (doca / 0.0625);

    double pdf_val = this->InterpolatePDF(time_residual, sigma, tau_eff);
    pdf_val *= doca_penalty; // unused

    if (pdf_val < 1e-3) {
      pdf_val = 1e-3;
    }
    llike -= log(pdf_val);
  }

  for (int i = 0; i < 7; i++) {
    if (this->constraints[i] > 0) {
      llike += pow((x[i] - this->constraint_means[i]) / this->constraints[i], 2);
    }
  }
  return (double)llike;
}

double GaussianDriftFit::operator()(const std::vector<double>& x) const {
  double const&a0 = x[0];
  double const&b0 = x[1];
  double const&a1 = x[2];
  double const&b1 = x[3];
  double t0 = x[4];
  long double llike = 0;

  CLHEP::Hep3Vector intercept(a0, 0, b0);
  CLHEP::Hep3Vector dir(a1, -1, b1);
  dir = dir.unit();

  std::vector<double> times(this->shs.size(),0);
  std::vector<double> tres(this->shs.size(),0);
  double average_time = 0;
  size_t count = 0;

  for (size_t i = 0; i < this->shs.size(); i++) {

    Straw const& straw = tracker->getStraw(this->shs[i].strawId());
    TwoLinePCA pca(intercept, dir, straw.getMidPoint(), straw.getDirection());

    if (constrainToStraw && pca.dca() > 2.5){
      llike += pow((pca.dca()-2.5)/0.1,2);
    }

    if (excludeHit == (int)i){
      continue;
    }
    count++;

    double longdist = (pca.point2() - straw.getMidPoint()).dot(straw.getDirection());
    if (fabs(longdist) > straw.halfLength()){
      llike += pow((fabs(longdist)-straw.halfLength()),2);
      longdist = std::copysign(straw.halfLength(),longdist);
    }
//    double longres = srep.wpRes(this->shs[i].energyDep() * 1000., fabs(longdist));
    double longres = srep.wpRes(this->shs[i].energyDep()*1000., this->shs[i].wireDist());

    llike += pow(longdist - this->shs[i].wireDist(), 2) / pow(longres, 2);

    double drift_time = srep.driftDistanceToTime(this->shs[i].strawId(), pca.dca(), 0) +
              srep.driftTimeOffset(this->shs[i].strawId(), pca.dca(), 0);

    double drift_res = srep.driftTimeError(this->shs[i].strawId(), pca.dca(), 0);

    double traj_time = ((pca.point1() - intercept).dot(dir)) / 299.9;
    double hit_t0 = this->shs[i].time() - this->shs[i].propTime() - traj_time - drift_time;
    //llike += pow((t0 - hit_t0) / drift_res, 2);
    average_time += hit_t0;
    times[i] = hit_t0;
    tres[i] = drift_res;

  }

  if (fixedT0){
    t0 = average_time/(int) count;
  }
  for (size_t i=0;i<this->shs.size();i++){
    if (excludeHit == (int)i){
      continue;
    }
    double drift_res = tres[i];
    if (fixedDriftRes)
      drift_res = driftRes;
    llike += pow((t0 - times[i]) / drift_res, 2);
  }
  //  std::cout << excludeHit << " " << fixedT0 << " " << t0 << " " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4] << "  =>  " <<
  //  llike << std::endl;

  return (double)llike;
}

double GaussianDriftFit::averageT0(const std::vector<double>& x) const {
  double const&a0 = x[0];
  double const&b0 = x[1];
  double const&a1 = x[2];
  double const&b1 = x[3];

  CLHEP::Hep3Vector intercept(a0, 0, b0);
  CLHEP::Hep3Vector dir(a1, -1, b1);
  dir = dir.unit();

  double average_time = 0;
  size_t count = 0;

  for (size_t i = 0; i < this->shs.size(); i++) {
    if (excludeHit == (int)i){
      continue;
    }
    count++;
    Straw const& straw = tracker->getStraw(this->shs[i].strawId());
    TwoLinePCA pca(intercept, dir, straw.getMidPoint(), straw.getDirection());

    double drift_time = srep.driftDistanceToTime(this->shs[i].strawId(), pca.dca(), 0) +
      srep.driftTimeOffset(this->shs[i].strawId(), pca.dca(), 0);

    double traj_time = ((pca.point1() - intercept).dot(dir)) / 299.9;
    double hit_t0 = this->shs[i].time() - this->shs[i].propTime() - traj_time - drift_time;
    average_time += hit_t0;
  }
  return average_time/(int) count;
}

double GaussianDriftFit::DOCAresidual(ComboHit const& sh, const std::vector<double>& x) const {
  double const& a0 = x[0];
  double const& b0 = x[1];
  double const& a1 = x[2];
  double const& b1 = x[3];
  double const& t0 = x[4];

  CLHEP::Hep3Vector intercept(a0, 0, b0);
  CLHEP::Hep3Vector dir(a1, -1, b1);
  dir = dir.unit();

  Straw const& straw = tracker->getStraw(sh.strawId());
  TwoLinePCA pca(intercept, dir, straw.getMidPoint(), straw.getDirection());
  double traj_time = ((pca.point1() - intercept).dot(dir)) / 299.9;

  double predictedDistance = pca.dca();
  double hit_t0 =
      sh.propTime() + traj_time + t0 + srep.driftTimeOffset(sh.strawId(), pca.dca(), 0);
  double measuredDistance = srep.driftTimeToDistance(sh.strawId(), sh.time() - hit_t0, 0);


  double resid = predictedDistance - measuredDistance;

  return resid;//(pca.s2() > 0 ? resid : -resid);
}

double GaussianDriftFit::reduced_chisq(const std::vector<double>& x) {
  double const& a0 = x[0];
  double const& b0 = x[1];
  double const& a1 = x[2];
  double const& b1 = x[3];

  CLHEP::Hep3Vector intercept(a0, 0, b0);
  CLHEP::Hep3Vector dir(a1, -1, b1);
  dir = dir.unit();

  double chi_sq = 0;
  size_t ndof = this->shs.size() - x.size();

  for (size_t i = 0; i < this->shs.size(); i++) {
    Straw const& straw = tracker->getStraw(shs[i].strawId());
    TwoLinePCA pca(intercept, dir, straw.getMidPoint(), straw.getDirection());

    double drift_res = srep.driftTimeError(shs[i].strawId(), pca.dca(), 0);
    chi_sq += pow(TimeResidual(shs[i], x) / drift_res, 2);
  }

  return chi_sq / ndof;
}

int GaussianDriftFit::HitAmbiguity(ComboHit const& sh, const std::vector<double>& x) const {
  double const& a0 = x[0];
  double const& b0 = x[1];
  double const& a1 = x[2];
  double const& b1 = x[3];

  CLHEP::Hep3Vector intercept(a0, 0, b0);
  CLHEP::Hep3Vector dir(a1, -1, b1);
  dir = dir.unit();

  Straw const& straw = tracker->getStraw(sh.strawId());

  Hep3Vector sep = intercept - straw.getMidPoint();
  Hep3Vector perp = (dir.cross(straw.getDirection())).unit();
  double dperp = perp.dot(sep);

  return (dperp > 0 ? -1 : 1);
}

double GaussianDriftFit::TimeResidual(ComboHit const& sh, const std::vector<double>& x) const {
  double const&a0 = x[0];
  double const&b0 = x[1];
  double const&a1 = x[2];
  double const&b1 = x[3];
  double t0 = x[4];

  CLHEP::Hep3Vector intercept(a0, 0, b0);
  CLHEP::Hep3Vector dir(a1, -1, b1);
  dir = dir.unit();

  Straw const& straw = tracker->getStraw(sh.strawId());
  TwoLinePCA pca(intercept, dir, straw.getMidPoint(), straw.getDirection());

  double drift_time = srep.driftDistanceToTime(sh.strawId(), pca.dca(), 0) +
    srep.driftTimeOffset(sh.strawId(), pca.dca(), 0);

  double traj_time = ((pca.point1() - intercept).dot(dir)) / 299.9;
  double hit_t0 = sh.time() - sh.propTime() - traj_time - drift_time;
  return t0 - hit_t0;
}

double GaussianDriftFit::DOCAresidualError(ComboHit const& sh, const std::vector<double>& x,
                                           const std::vector<double>& cov) const {
  double const& a0 = x[0];
  double const& b0 = x[1];
  double const& a1 = x[2];
  double const& b1 = x[3];
  //  double t0 = x[4];

  CLHEP::Hep3Vector intercept(a0, 0, b0);
  CLHEP::Hep3Vector dir(a1, -1, b1);
  dir = dir.unit();

  Straw const& straw = tracker->getStraw(sh.strawId());
  TwoLinePCA pca(straw.getMidPoint(), straw.getDirection(), intercept, dir);

  auto s0 = straw.getMidPoint();
  auto s1 = straw.getDirection();
  double x0 = s0.x();
  double y0 = s0.y();
  double z0 = s0.z();
  double x1 = s1.x();
  double y1 = s1.y();
  double z1 = s1.z();

  double da0 = (-(b1 * y1) + z1) /
               sqrt(pow(x1 + a1 * y1, 2) + pow(-b1 * y1 - z1, 2) + pow(b1 * x1 - a1 * z1, 2));
  double db0 = (a1 * y1 + x1) /
               sqrt(pow(x1 + a1 * y1, 2) + pow(-b1 * y1 - z1, 2) + pow(b1 * x1 - a1 * z1, 2));
  double da1 =
      -(((x1 + a1 * y1) * (b0 - z0) + (a0 - x0) * (-(b1 * y1) + z1) - y0 * (b1 * x1 - a1 * z1)) *
        (2 * y1 * (x1 + a1 * y1) - 2 * z1 * (b1 * x1 - a1 * z1))) /
          (2 * pow(pow(x1 + a1 * y1, 2) + pow(-(b1 * y1) - z1, 2) + pow(b1 * x1 - a1 * z1, 2),
                   3 / 2.)) +
      (y1 * (b0 - z0) + y0 * z1) /
          sqrt(pow(x1 + a1 * y1, 2) + pow(-(b1 * y1) - z1, 2) + pow(b1 * x1 - a1 * z1, 2));
  double db1 =
      -((-2 * y1 * (-(b1 * y1) - z1) + 2 * x1 * (b1 * x1 - a1 * z1)) *
        ((x1 + a1 * y1) * (b0 - z0) + (a0 - x0) * (-(b1 * y1) + z1) - y0 * (b1 * x1 - a1 * z1))) /
          (2 * pow(pow(x1 + a1 * y1, 2) + pow(-(b1 * y1) - z1, 2) + pow(b1 * x1 - a1 * z1, 2),
                   3 / 2.)) +
      (-(x1 * y0) - (a0 - x0) * y1) /
          sqrt(pow(x1 + a1 * y1, 2) + pow(-(b1 * y1) - z1, 2) + pow(b1 * x1 - a1 * z1, 2));
  double dt0 = srep.driftInstantSpeed(sh.strawId(), pca.dca(), 0);

  std::vector<double> dx = {da0, db0, da1, db1, dt0};
  std::vector<double> temp(5, 0);
  for (size_t i = 0; i < 5; i++) {
    for (size_t j = 0; j < 5; j++) {
      double thiscov = cov[i + j * (j + 1) / 2];
      if (j < i)
        thiscov = cov[j + i * (i + 1) / 2];

      temp[i] += dx[j] * thiscov;
    }
  }

  double residerror = 0;
  for (size_t i = 0; i < 5; i++) {
    residerror += temp[i] * dx[i];
  }

  return sqrt(residerror);
}
