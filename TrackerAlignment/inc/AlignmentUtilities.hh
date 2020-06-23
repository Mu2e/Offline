// Ryunosuke O'Neil, 2020
// roneil@fnal.gov
// ryunoneil@gmail.com

// ROOT
#include "CLHEP/Vector/ThreeVector.h"
#include "DbTables/inc/TrkAlignPlane.hh"
#include "GeneralUtilities/inc/HepTransform.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "RtypesCore.h"
#include "TTree.h"

// art
#include "TrackerConditions/inc/StrawResponse.hh"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

// Offline

#include "DbTables/inc/TrkAlignPanel.hh"

#include "DataProducts/inc/StrawId.hh"

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"

#include "TrackerAlignment/inc/AlignmentDerivatives.hh"
#include "TrackerAlignment/inc/MilleDataWriter.hh"
#include <cstdint>
#include <iterator>
#include <vector>

using namespace mu2e;

namespace AlignmentUtilities {


// KinKal-inspired track description
struct CosmicTimeTrack {
  double params[5];

  enum paraml { a0 = 0, b0 = 1, a1 = 2, b1 = 3, t0 = 4 };

  CosmicTimeTrack(double a0, double b0, double a1, double b1, double t0) :
      params{a0, b0, a1, b1, t0} {}
  
  CosmicTimeTrack(std::vector<double> const& pars) :
      params{pars[a0],pars[b0], pars[a1], pars[b1], pars[t0]} {}

  std::vector<double> as_vector() {
    return {params[a0], params[b0], params[a1], params[b1], params[t0]};
  }

  Hep3Vector intercept() const { 
    return {params[a0], 0, params[b0]}; 
  }

  Hep3Vector direction() const {
    Hep3Vector result{params[a1], -1, params[b1]};
    return result.unit();
  }

  size_t npars() const { return 5; }

  void setParams(std::vector<double> const& pars) {
    params[a0] = pars[a0];
    params[a1] = pars[a1];
    params[b0] = pars[b0];
    params[b1] = pars[b1];
    params[t0] = pars[t0];
  }
  void setParams(CosmicTimeTrack const& pars) {
    params[a0] = pars.params[a0];
    params[a1] = pars.params[a1];
    params[b0] = pars.params[b0];
    params[b1] = pars.params[b1];
    params[t0] = pars.params[t0];
  }
};

bool testDerivatives(
    TwoLinePCA const& expected_pca,
    Tracker const& alignedTracker,

    CosmicTimeTrack const& track,
    StrawId const& strawId,
    TrkAlignPlane::Row const&rowpl,
    TrkAlignPanel::Row const&rowpa,
    Tracker const& nominalTracker,
    StrawResponse const& strawRes);


std::pair<std::vector<double>, std::vector<double>> 
  analyticalDerivatives(CosmicTimeTrack const& track,
    StrawId const& strawId,
    TrkAlignPlane::Row const&rowpl,
    TrkAlignPanel::Row const&rowpa,
    Tracker const& nominalTracker,
    double const& driftvel);

TMatrixD residualCovariance(CosmicTimeTrack const& track, 
  std::vector<double> const& track_cov,
  std::vector<std::vector<double>> const& local_derivatives,
  std::vector<double> const& meas_err);

int hitAmbiguity(CosmicTimeTrack const& track, Hep3Vector const& straw_mp,
                 Hep3Vector const& straw_dir);

/* diagnostic prints
 */

void diagPrintTrack(CosmicTimeTrack const& track);

void diagPrintHit(CosmicTimeTrack const& track,
  double const& tresid, double const& resolution,
  std::vector<double> const& derivativesLocal,
  std::vector<double> const& derivativesGlobal,
  StrawId const& strawId);


std::pair<Hep3Vector, Hep3Vector> alignStraw(Tracker const& tracker, Plane const& plane,
                                             Panel const& panel, StrawId const& strawId,
                                             HepTransform const& align_tracker,
                                             HepTransform const& align_plane,
                                             HepTransform const& align_panel);

/* Numerical partial DOCA/TOCA derivatives
 *
 */ 

std::pair<std::vector<double>, std::vector<double>>
numericalDerivatives(CosmicTimeTrack const& _track, StrawId const& straw,
                         TrkAlignPlane::Row const& alignPlane,
                         TrkAlignPanel::Row const& alignPanel,
                         Tracker const& nominalTracker, 
                         StrawResponse const& strawRes, 
                         bool useTimeDomain = true);

} // namespace AlignmentUtilities
