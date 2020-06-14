// Ryunosuke O'Neil, 2019
// roneil@fnal.gov
// ryunoneil@gmail.com

// An helper class for carrying out Alignment diagnostics

// ROOT
#include "CLHEP/Vector/ThreeVector.h"
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

#include "TrackerAlignment/inc/AlignStraw.hh"
#include "TrackerAlignment/inc/AlignmentDerivatives.hh"
#include "TrackerAlignment/inc/MilleDataWriter.hh"
#include <cstdint>
#include <iterator>
#include <vector>

using namespace mu2e;

namespace AlignmentUtilities {

struct CosmicTimeTrackData {
    double 
}



double calcDCA(AlignTrackParams const& trackParams, StrawId const& strawId,
                                      std::vector<double> const& globals) const {
  double const& a0 = trackParams[0];
  double const& b0 = trackParams[1];
  double const& a1 = trackParams[2];
  double const& b1 = trackParams[3];
  // double const& t0 = trackParams[4]; // unused

  CLHEP::Hep3Vector intercept(a0, 0, b0);
  CLHEP::Hep3Vector dir(a1, -1, b1);
  dir = dir.unit();

  Plane const& nominal_plane = nominalTracker.getPlane(strawId);
  Panel const& nominal_panel = nominalTracker.getPanel(strawId);

  HepTransform align_tracker{0, 0, 0, 0, 0, 0};

  HepTransform align_plane{globals[0], globals[1], globals[2], globals[3], globals[4], globals[5]};

  HepTransform align_panel{globals[6], globals[7],  globals[8],
                           globals[9], globals[10], globals[11]};

  // returns pair of vectors { straw_pos, straw_dir }
  auto aligned_result = AlignStraw::alignStraw(nominalTracker, nominal_plane, nominal_panel,
                                               strawId, align_tracker, align_plane, align_panel);

  TwoLinePCA pca(intercept, dir, aligned_result.first, aligned_result.second);

  double result = (pca.s2() > 0 ? pca.dca() : -pca.dca());

  return result;
}

double AlignmentUtilities::numericalDerivative(StrawId const& straw, AlignTrackParams& local,
                                                  std::vector<double>& global, bool isGlobalParam,
                                                  int const& paramIdx, double step_size) const {
  // calculate partial derivative wrt param at paramIdx in either
  // local, or global param array

  if (isGlobalParam) {
    global[paramIdx] -= step_size;
  } else {
    local[paramIdx] -= step_size;
  }

  double pdiff = strawRes.driftDistanceToTime(straw, calcDCA(local, straw, global), 0);

  if (isGlobalParam) {
    global[paramIdx] += 2.0 * step_size;
  } else {
    local[paramIdx] += 2.0 * step_size;
  }

  pdiff -= strawRes.driftDistanceToTime(straw, calcDCA(local, straw, global), 0);
  pdiff /= (2.0 * step_size);

  return pdiff;
}

std::pair<std::vector<double>, std::vector<double>>
AlignTrackDiagnostics::calcNumericalDerivatives(StrawId const& straw) const {
  std::vector<double> result_locals;
  std::vector<double> result_globals;

  // construct for current track, and track hit

  // we copy here because we modify the contents as we calculate each derivative
  AlignTrackParams track(currentTrackParams);

  auto const& alignPlane = alignmentConstants_plane[straw.getPlane()];
  auto const& alignPanel = alignmentConstants_panel[straw.getPanel()];

  std::vector<double> globals(alignPlane);
  globals.insert(globals.end(), alignPanel.begin(), alignPanel.end());

  // locals first
  for (size_t paramIdx = 0; paramIdx < track.size(); ++paramIdx) {
    result_locals.emplace_back(numericalDerivative(straw, track, globals, false, paramIdx, 1e-7));
  }

  for (size_t paramIdx = 0; paramIdx < globals.size(); ++paramIdx) {
    result_globals.emplace_back(numericalDerivative(straw, track, globals, true, paramIdx, 1e-7));
  }

  std::pair<std::vector<double>, std::vector<double>> result{result_locals, result_globals};

  return result;
}

} // namespace AlignmentUtilities