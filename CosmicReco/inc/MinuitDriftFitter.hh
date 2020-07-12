#ifndef _COSMIC_RECO_MINUITDRIFTFITTER_HH
#define _COSMIC_RECO_MINUITDRIFTFITTER_HH

#include "CosmicReco/inc/PDFFit.hh"
#include "DataProducts/inc/XYZVec.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "TrackerConditions/inc/StrawDrift.hh"
#include "TrackerGeom/inc/Tracker.hh"

// For Drift:
#include "BTrk/BaBar/BaBar.hh"
#include "BTrk/BbrGeom/BbrVectorErr.hh"
#include "BTrk/BbrGeom/HepPoint.h"
#include "BTrk/BbrGeom/Trajectory.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "BTrk/ProbTools/ChisqConsistency.hh"
#include "BTrk/TrkBase/TrkMomCalculator.hh"
#include "BTrk/TrkBase/TrkPoca.hh"
#include "BTrkData/inc/TrkStrawHit.hh"

// ROOT
#include "TF1.h"
#include "TH1F.h"
#include "TMath.h"

// Minuit
#include <Minuit2/FCNBase.h>

using namespace mu2e;
struct FitResult {
public:
  std::vector<std::string> names;
  std::vector<double> bestfit;
  std::vector<double> bestfiterrors;
  std::vector<double> bestfitcov;

  std::vector<double> StartDOCAs;
  std::vector<double> StartTimeResiduals;

  std::vector<double> GaussianEndDOCAs;
  std::vector<double> GaussianEndTimeResiduals;

  std::vector<double> FullFitEndDOCAs;
  std::vector<double> FullFitEndTimeResiduals;

  std::vector<double> RecoAmbigs;

  double NLL;
};

namespace MinuitDriftFitter {
FitResult DoFit(int const& _diag, CosmicTrackSeed& tseed, StrawResponse const& srep,
                const Tracker* tracker, double const& max_doca, unsigned int const& minChits,
                int const& MaxLogL, double const& _gaussTres, double const& maxTres);


void DoDriftTimeFit(
    std::vector<double> & pars, 
    std::vector<double> & errors,
    std::vector<double> & cov_out,
    bool & minuit_converged,
    GaussianDriftFit const& fit, 
    int diag=0, double mntolerance=0.1, double mnprecision=-1);

void DoDriftTimeFit(int const& diag, CosmicTrackSeed& tseed, StrawResponse const& srep,
                    const Tracker* tracker, double mntolerance=0.1, double mnprecision=-1);

} // namespace MinuitDriftFitter

#endif
