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

// could change in future
typedef std::vector<double> AlignTrackParams;

class AlignTrackDiagnostics {

private:

  static constexpr int MAX_NHITS = 100;

  // handles passed from art module
  art::ServiceHandle<art::TFileService> const& fileService;
  Tracker const& nominalTracker;
  Tracker const& alignedTracker;
  StrawResponse const& strawRes;

  std::vector<std::vector<double>> alignmentConstants_plane;
  std::vector<std::vector<double>> alignmentConstants_panel;

  AlignTrackParams currentTrackParams;

  // Alignment constant TTree
  TTree* alignconst_tree;

  Int_t object_type; // 1 = plane, 2 = panel
  Int_t object_id;   // plane 0-35 or panel 0,215, etc

  // shifts
  Float_t dx;
  Float_t dy;
  Float_t dz;

  // rotations
  Float_t rx;
  Float_t ry;
  Float_t rz;

  // Track diagnostic TTree
  TTree* trackdiag_tree;

  Int_t nHits;
  Float_t doca_residual[MAX_NHITS];
  Float_t time_residual[MAX_NHITS];
  Float_t residual_err[MAX_NHITS];

  Float_t drift_res_doca[MAX_NHITS];
  Float_t drift_res_time[MAX_NHITS];

  Float_t doca_pulls[MAX_NHITS];
  Float_t time_pulls[MAX_NHITS];

  Float_t doca[MAX_NHITS];
  Float_t time[MAX_NHITS];

  Int_t plane_id[MAX_NHITS];
  Int_t panel_id[MAX_NHITS];

  Double_t A0;
  Double_t A1;
  Double_t B0;
  Double_t B1;
  Double_t T0;

  Double_t chisq;
  Double_t chisq_doca;
  Int_t ndof;
  Double_t pvalue;

  Int_t panels_trav;
  Int_t planes_trav;

  Int_t track_ambiguity_s2;
  Int_t track_ambiguity_utilfn;

  void initTrackTree();
  void initAlignConstTree();

  // other

  double numericalDerivative(StrawId const& straw,
    AlignTrackParams & local,
    std::vector<double> & global,
    bool isGlobalParam,
    int const& paramIdx, 
    double step_size = 1e-7) const;

  std::pair<std::vector<double>, std::vector<double>> 
    calcNumericalDerivatives(
      StrawId const& straw) const;
      
  double calcDCA(
    AlignTrackParams const& trackParams,
    StrawId const& strawId,
    std::vector<double> const& globals) const;

public:
  AlignTrackDiagnostics(
    art::ServiceHandle<art::TFileService> const& tfs, 
    Tracker const& nomTracker, 
    Tracker const& alignTracker,
    StrawResponse const& strawRes) : 
      fileService(tfs), 
      nominalTracker(nomTracker), 
      alignedTracker(alignTracker),
      strawRes(strawRes) {}

  // write alignment constants used
  void writeConstants(TrkAlignPanel const& alignConstPlanes, TrkAlignPanel const& alignConstPanels);

  // call when a new track is found
  void newTrack(AlignTrackParams const& params);

  // call when finished processing a track
  void finishTrack(bool isGoodTrack);

  // call when a track hit is processed
  // (after newTrack is called, and before finishTrack is called)
  void newTrackHit();

  void printTrack();
};

void AlignTrackDiagnostics::initTrackTree() {
  trackdiag_tree =
      fileService->make<TTree>("tracks", "Tracks collected for an alignment iteration");
  trackdiag_tree->Branch("nHits", &nHits, "nHits/I");
  trackdiag_tree->Branch("doca_resid", &doca_residual, "doca_resid[nHits]/F");
  trackdiag_tree->Branch("time_resid", &time_residual, "time_resid[nHits]/F");
  trackdiag_tree->Branch("drift_res_doca", &drift_res_doca, "drift_res_doca[nHits]/F");
  trackdiag_tree->Branch("drift_res", &drift_res_time, "drift_res[nHits]/F");

  trackdiag_tree->Branch("resid_err", &residual_err, "resid_err[nHits]/F");

  trackdiag_tree->Branch("pull_doca", &doca_pulls, "pull_doca[nHits]/F");
  trackdiag_tree->Branch("pull_hittime", &time_pulls, "pull_doca[nHits]/F");

  trackdiag_tree->Branch("doca", &doca, "doca[nHits]/F");
  trackdiag_tree->Branch("time", &time, "time[nHits]/F");
  trackdiag_tree->Branch("plane", &plane_id, "plane[nHits]/I");
  trackdiag_tree->Branch("panel", &panel_id, "panel[nHits]/I");

  trackdiag_tree->Branch("A0", &A0, "A0/D");
  trackdiag_tree->Branch("A1", &A1, "A1/D");
  trackdiag_tree->Branch("B0", &B0, "B0/D");
  trackdiag_tree->Branch("B1", &B1, "B1/D");
  trackdiag_tree->Branch("T0", &T0, "T0/D");

  trackdiag_tree->Branch("chisq", &chisq, "chisq/D");
  trackdiag_tree->Branch("chisq_doca", &chisq_doca, "chisq_doca/D");

  trackdiag_tree->Branch("ndof", &ndof, "ndof/I");
  trackdiag_tree->Branch("pvalue", &pvalue, "pvalue/D");

  trackdiag_tree->Branch("panels_trav", &panels_trav, "panels_trav/I");
  trackdiag_tree->Branch("planes_trav", &planes_trav, "planes_trav/I");
}

void AlignTrackDiagnostics::initAlignConstTree() {
  alignconst_tree =
      fileService->make<TTree>("alignconst", "Alignment constants for this iteration");
  alignconst_tree->Branch("object_type", &nHits, "object_type/I");
  alignconst_tree->Branch("object_id", &doca_residual, "object_id/I");

  alignconst_tree->Branch("dx", &time_residual, "dx/F");
  alignconst_tree->Branch("dy", &drift_res_doca, "dy/F");
  alignconst_tree->Branch("dz", &drift_res_time, "dz/F");

  alignconst_tree->Branch("rx", &time_residual, "rx/F");
  alignconst_tree->Branch("ry", &drift_res_doca, "ry/F");
  alignconst_tree->Branch("rz", &drift_res_time, "rz/F");
}

double AlignTrackDiagnostics::calcDCA(
  AlignTrackParams const& trackParams,
  StrawId const& strawId,
  std::vector<double> const& globals) const
{
  double const& a0 = trackParams[0];
  double const& b0 = trackParams[1];
  double const& a1 = trackParams[2];
  double const& b1 = trackParams[3];
  //double const& t0 = trackParams[4]; // unused

  CLHEP::Hep3Vector intercept(a0, 0, b0);
  CLHEP::Hep3Vector dir(a1, -1, b1);
  dir = dir.unit();

  Plane const& nominal_plane = nominalTracker.getPlane(strawId);
  Panel const& nominal_panel = nominalTracker.getPanel(strawId);

  HepTransform align_tracker { 0,0,0, 0,0,0 };

  HepTransform align_plane { 
    globals[0], globals[1], globals[2],
    globals[3], globals[4], globals[5] };

  HepTransform align_panel { 
    globals[6], globals[7], globals[8],
    globals[9], globals[10], globals[11] };

  // returns pair of vectors { straw_pos, straw_dir }
  auto aligned_result = AlignStraw::alignStraw(nominalTracker, nominal_plane, 
      nominal_panel, strawId, align_tracker, align_plane, align_panel);

  TwoLinePCA pca(intercept, dir, aligned_result.first, aligned_result.second);

  double result = (pca.s2() > 0 ? pca.dca() : -pca.dca());

  return result;
}


double AlignTrackDiagnostics::numericalDerivative(StrawId const& straw,
                                                AlignTrackParams & local,
                                                std::vector<double> & global,
                                                bool isGlobalParam,
                                                int const& paramIdx, 
                                                double step_size) const {
  // calculate partial derivative wrt param at paramIdx in either 
  // local, or global param array

  if (isGlobalParam) {
    global[paramIdx] -= step_size;
  }
  else {
    local[paramIdx] -= step_size;
  }

  double pdiff = strawRes.driftDistanceToTime(
            straw,
            calcDCA(local, straw, global),
            0);

  if (isGlobalParam) {
    global[paramIdx] += 2.0*step_size;
  }
  else {
    local[paramIdx] += 2.0*step_size;
  }
    
  pdiff -= strawRes.driftDistanceToTime(
              straw,
              calcDCA(local, straw, global),
              0);
  pdiff /= (2.0 * step_size);

  return pdiff;
}

std::pair<std::vector<double>, std::vector<double>> 
  AlignTrackDiagnostics::calcNumericalDerivatives(
    StrawId const& straw) const {
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
    result_locals.emplace_back(
      numericalDerivative(straw,
        track,
        globals,
        false,
        paramIdx,
        1e-7)
    );
  }

  for (size_t paramIdx = 0; paramIdx < globals.size(); ++paramIdx) {
    result_globals.emplace_back(
      numericalDerivative(straw,
        track,
        globals,
        true,
        paramIdx,
        1e-7)
    );
  }

  std::pair<std::vector<double>, std::vector<double>> result {
    result_locals,
    result_globals
  };


  return result;
}