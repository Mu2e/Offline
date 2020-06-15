// Ryunosuke O'Neil, 2020
// roneil@fnal.gov
// ryunoneil@gmail.com

#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>
#include "RtypesCore.h"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "DbTables/inc/TrkAlignPlane.hh"
#include "GeneralUtilities/inc/HepTransform.hh"
#include "Minuit2/MnUserCovariance.h"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "TTree.h"
#include "DbTables/inc/TrkAlignPanel.hh"
#include "DataProducts/inc/StrawId.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerAlignment/inc/AlignmentDerivatives.hh"
#include "TrackerAlignment/inc/MilleDataWriter.hh"
#include "TrackerAlignment/inc/AlignmentUtilities.hh"

using namespace mu2e;

namespace AlignmentUtilities {

bool testDerivatives(
    TwoLinePCA const& expected_pca,
    Tracker const& alignedTracker,

    CosmicTimeTrack const& track,
    StrawId const& strawId,
    TrkAlignPlane::Row const&rowpl,
    TrkAlignPanel::Row const&rowpa,
    Tracker const& nominalTracker,
    StrawResponse const& strawRes) {
  double tolerance = 1e-7;

  // first compare DOCA
  Straw const& alignedStraw = alignedTracker.getStraw(strawId);
  double expected_dca = hitAmbiguity(track, 
      alignedStraw.getMidPoint(), 
      alignedStraw.getDirection()) * expected_pca.dca();

  Straw const& nominalStraw = nominalTracker.getStraw(strawId);
  Hep3Vector const& nominalStraw_mp = nominalStraw.getMidPoint();
  Hep3Vector const& nominalStraw_dir = nominalStraw.getDirection();
  Hep3Vector const& plane_origin = nominalTracker.getPlane(strawId).origin();
  Hep3Vector const& panel_origin = nominalTracker.getPanel(strawId).straw0MidPoint();
  
  // calculate DOCA from analytical generated code
  double doca = CosmicTrack_DCA(
      track.params[CosmicTimeTrack::a0], 
      track.params[CosmicTimeTrack::b0], 
      track.params[CosmicTimeTrack::a1], 
      track.params[CosmicTimeTrack::b1], 
      track.params[CosmicTimeTrack::t0], 
      rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(), rowpl.ry(), rowpl.rz(), 
      rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(), rowpa.ry(), rowpa.rz(),

      nominalStraw_mp.x(), nominalStraw_mp.y(), nominalStraw_mp.z(), nominalStraw_dir.x(),
      nominalStraw_dir.y(), nominalStraw_dir.z(), plane_origin.x(), plane_origin.y(),
      plane_origin.z(), panel_origin.x(), panel_origin.y(), panel_origin.z(), 0);

  if (std::abs(doca - expected_dca) > tolerance) {
    std::cerr << "doca mismatch: diff = " << std::abs(doca - expected_dca) << std::endl;
    return false;
  }

  double driftvel = strawRes.driftInstantSpeed(strawId, doca, 0);
  auto anaDerivatives = analyticalDerivatives(track, strawId, rowpl, rowpa, nominalTracker, driftvel);
  auto numDerivatives = numericalDerivatives(track, strawId, rowpl, rowpa, nominalTracker, strawRes);

  // compare local derivatives

  if (numDerivatives.first.size() != anaDerivatives.first.size()) {
    std::cerr << "size mismatch  (local)" << std::endl;
    return false;
  }
  
  for (size_t i = 0; i < numDerivatives.first.size(); ++i) {
     if (std::abs(anaDerivatives.first[i] - numDerivatives.first[i]) > tolerance) {
      std::cerr << "local derivative mismatch: diff = " 
                << std::abs(anaDerivatives.first[i] - numDerivatives.first[i]) 
                << std::endl;
       return false;
     }
  }

  // compare global derivatives
  if (numDerivatives.second.size() != anaDerivatives.second.size()) {
    std::cerr << "size mismatch  (global)" << std::endl;
    return false;
  }
  for (size_t i = 0; i < numDerivatives.second.size(); ++i) {
     if (std::abs(anaDerivatives.second[i] - numDerivatives.second[i]) < tolerance) {
       std::cerr << "global derivative mismatch: diff = " 
                << std::abs(anaDerivatives.second[i] - numDerivatives.second[i]) 
                << std::endl;
       return false;
     }
  }

  return true;
}

std::pair<std::vector<double>, std::vector<double>> 
  analyticalDerivatives(CosmicTimeTrack const& track,
    StrawId const& strawId,
    TrkAlignPlane::Row const&rowpl,
    TrkAlignPanel::Row const&rowpa,
    Tracker const& nominalTracker,
    double const& driftvel) {

  auto const& plane_origin = nominalTracker.getPlane(strawId.getPlane()).origin();
  auto const& panel_origin = nominalTracker.getPanel(strawId).straw0MidPoint();

  Straw const& nominalStraw = nominalTracker.getStraw(strawId);
  auto const& nominalStraw_mp = nominalStraw.getMidPoint();
  auto const& nominalStraw_dir = nominalStraw.getDirection();

  auto derivativesLocal = CosmicTrack_DCA_LocalDeriv(
      track.params[CosmicTimeTrack::a0], 
      track.params[CosmicTimeTrack::b0], 
      track.params[CosmicTimeTrack::a1], 
      track.params[CosmicTimeTrack::b1], 
      track.params[CosmicTimeTrack::t0], 
      rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(), rowpl.ry(), rowpl.rz(), 
      rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(), rowpa.ry(), rowpa.rz(),

      nominalStraw_mp.x(), nominalStraw_mp.y(), nominalStraw_mp.z(), nominalStraw_dir.x(),
      nominalStraw_dir.y(), nominalStraw_dir.z(), plane_origin.x(), plane_origin.y(),
      plane_origin.z(), panel_origin.x(), panel_origin.y(), panel_origin.z(), driftvel);
  
  auto derivativesGlobal = CosmicTrack_DCA_GlobalDeriv(
      track.params[CosmicTimeTrack::a0], 
      track.params[CosmicTimeTrack::b0], 
      track.params[CosmicTimeTrack::a1], 
      track.params[CosmicTimeTrack::b1], 
      track.params[CosmicTimeTrack::t0],  
  
      rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(), rowpl.ry(), rowpl.rz(), 
      rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(), rowpa.ry(), rowpa.rz(),

      nominalStraw_mp.x(), nominalStraw_mp.y(), nominalStraw_mp.z(), 
      nominalStraw_dir.x(),nominalStraw_dir.y(), nominalStraw_dir.z(),
      plane_origin.x(), plane_origin.y(), plane_origin.z(), 
      panel_origin.x(), panel_origin.y(), panel_origin.z(), 
      driftvel);

  return {derivativesLocal, derivativesGlobal};
}

int hitAmbiguity(CosmicTimeTrack const& track, Hep3Vector const& straw_mp,
                 Hep3Vector const& straw_dir) {
  Hep3Vector sep = track.intercept() - straw_mp;
  Hep3Vector perp = (track.direction().cross(straw_dir)).unit();
  double dperp = perp.dot(sep);

  return (dperp > 0 ? -1 : 1);
}

// cov(r) = V - HCH
TMatrixD residualCovariance(CosmicTimeTrack const& track, 
  std::vector<double> const& track_cov,
  std::vector<std::vector<double>> const& local_derivatives,
  std::vector<double> const& meas_err) {
  size_t const& nHits = meas_err.size();

  // V = cov( measurements )
  TMatrixD meas_cov(nHits, nHits);
  {
    meas_cov.Zero();
    for (size_t i = 0; i < nHits; i++) {
        meas_cov(i, i) = meas_err[i]*meas_err[i];
    }
  }

  // H = partial d (residuals)/ d (track params)
  TMatrixD residDerivs(nHits, 5);
  {
    residDerivs.Zero();
    for (size_t r = 0; r < nHits; r++) {
      for (size_t c = 0; c < 5; c++) {
        residDerivs(r, c) = local_derivatives[r][c];
      }
    }
  }

  // C = cov ( track parameters )
  TMatrixD trackcov_mat(5, 5);
  {
    ROOT::Minuit2::MnUserCovariance cov(track_cov, 5);
    for (size_t r = 0; r < 5; r++) {
      for (size_t c = 0; c < 5; c++) {
        trackcov_mat(r, c) = cov(r, c);
      }
    }
  }


  // meas_cov = residual_cov = V - H C H^T
  {
    TMatrixD HCH(nHits, nHits);
    TMatrixD HC(nHits, 5);
    HC.Mult(residDerivs, trackcov_mat);

    // Using MultT so HCH^T = HC * H^T with H = resid_local_derivs
    HCH.MultT(HC, residDerivs);
    meas_cov -= HCH;
  }

  return meas_cov;
}

/* 
 * diagnostic prints
 */

void diagPrintTrack(CosmicTimeTrack const& track) {
  std::cout << "[A0,B0,A1,B1,T0]: [" 
        <<         track.params[track.a0] 
        << ", " << track.params[track.b0] 
        << ", " << track.params[track.a1] 
        << ", " << track.params[track.b1]  
        << ", " << track.params[track.t0];
  std::cout << "]" << std::endl;
}

void diagPrintHit(CosmicTimeTrack const& track,
  double const& tresid, double const& resolution,
  std::vector<double> const& derivativesLocal,
  std::vector<double> const& derivativesGlobal,
  StrawId const& strawId) {
  std::cout << "TOCA resid: " 
            << tresid 
            << " +- " 
            << resolution 
            << std::endl;
  std::cout << "dr/d([A0,B0,A1,B1,T0]): [" 
            <<         derivativesLocal[0] 
            << ", " << derivativesLocal[1]                    
            << ", " << derivativesLocal[2]                     
            << ", " << derivativesLocal[3]                 
            << ", " << derivativesLocal[4];               
  std::cout << "]" << std::endl;

  std::cout << "dr/d(plane" 
            << strawId.getPlane() 
            << "[x, y, z]): [" 
            <<         derivativesGlobal[0] 
            << ", " << derivativesGlobal[1] 
            << ", " << derivativesGlobal[2];
  std::cout << "]" << std::endl;

}

// straw alignment function
// carbon copy of AlignedTrackerMaker::fromDb
std::pair<Hep3Vector, Hep3Vector> alignStraw(Tracker const& tracker, Plane const& plane,
                                             Panel const& panel, StrawId const& strawId,
                                             HepTransform const& align_tracker,
                                             HepTransform const& align_plane,
                                             HepTransform const& align_panel) {
  std::pair<Hep3Vector, Hep3Vector> result;
  
  // the whole tracker has nominal center on 0,0,0
  // how to place the plane in the tracker
  HepTransform plane_to_tracker(0.0, 0.0, plane.origin().z(), 0.0, 0.0, 0.0);

  // make an intermediate multiplication
  HepTransform plane_temp = align_tracker * (plane_to_tracker * align_plane);

  // how to place the panel in the plane
  Hep3Vector dv = panel.straw0MidPoint() - plane_to_tracker.displacement();
  double rz = dv.phi();

  HepTransform panel_to_plane(dv.x(), dv.y(), dv.z(), 0.0, 0.0, rz);

  // make an intermediate multiplication
  HepTransform panel_temp = plane_temp * (panel_to_plane * align_panel);

  Straw const& straw = tracker.getStraw(strawId);

  // how to place the straw in the panel
  double dx = straw.getMidPoint().perp() - panel.straw0MidPoint().perp();
  double dz = (straw.getMidPoint() - panel.straw0MidPoint()).z();

  Hep3Vector straw_to_panel = Hep3Vector(dx, 0.0, dz);
  Hep3Vector straw_dir = Hep3Vector(0.0, 1.0, 0.0);

  Hep3Vector aligned_straw = panel_temp * straw_to_panel;
  Hep3Vector aligned_straw_dir = panel_temp.rotation() * straw_dir;

  // aligned straw position inserted in the Tracker object

  Hep3Vector pdif = aligned_straw - straw.getMidPoint();
  Hep3Vector ddif = aligned_straw_dir - straw.getDirection();

  result.first = aligned_straw;
  result.second = aligned_straw_dir;

  return result;
}



/* Numerical partial DOCA/TOCA derivatives
 *
 */ 


namespace {

double docaGlobalDep(CosmicTimeTrack const& track, StrawId const& strawId,
                              std::vector<double> const& globals, Tracker const& nominalTracker) {

  Plane const& nominal_plane = nominalTracker.getPlane(strawId);
  Panel const& nominal_panel = nominalTracker.getPanel(strawId);

  HepTransform align_tracker{0, 0, 0, 0, 0, 0};
  HepTransform align_plane{globals[0], globals[1], globals[2], globals[3], globals[4], globals[5]};

  HepTransform align_panel{globals[6], globals[7],  globals[8],
                           globals[9], globals[10], globals[11]};

  // returns pair of vectors { straw_pos, straw_dir }
  auto aligned_result = alignStraw(nominalTracker, nominal_plane, nominal_panel,
                                               strawId, align_tracker, align_plane, align_panel);

  TwoLinePCA pca(track.intercept(), track.direction(), aligned_result.first, aligned_result.second);

  int ambig = hitAmbiguity(track, aligned_result.first, aligned_result.second);
  double result = ambig * pca.dca();

  return result;
}

// not meant to be called from outside of this namespace
double _numericalDerivative(StrawId const& straw, CosmicTimeTrack& track,
                           std::vector<double>& globals, Tracker const& nominalTracker,
                           StrawResponse const& strawRes, bool isGlobalParam,
                           size_t const& paramIdx, double step_size) {
  // calculate numerical partial derivative wrt param at paramIdx in either
  // local, or global param array

  if (isGlobalParam) {
    globals[paramIdx] -= step_size;
  } else {
    track.params[paramIdx] -= step_size;
  }

  // double pdiff = strawRes.driftDistanceToTime(
  //     straw, docaGlobalDep(track, straw, globals, nominalTracker), 0);
  double pdiff = docaGlobalDep(track, straw, globals, nominalTracker);
  double driftvel = strawRes.driftInstantSpeed(straw, pdiff, 0);
  pdiff /= driftvel;

  if (isGlobalParam) {
    globals[paramIdx] += 2.0 * step_size;
  } else {
    track.params[paramIdx] += 2.0 * step_size;
  }

  // pdiff -= strawRes.driftDistanceToTime(
  //     straw, docaGlobalDep(track, straw, globals, nominalTracker), 0);
  double doca2 = docaGlobalDep(track, straw, globals, nominalTracker);
  driftvel = strawRes.driftInstantSpeed(straw, pdiff, 0);
  pdiff -= doca2 / driftvel;

  pdiff /= (2.0 * step_size);

  return pdiff;
}
}

std::pair<std::vector<double>, std::vector<double>>
numericalDerivatives(CosmicTimeTrack const& _track, StrawId const& straw,
                         TrkAlignPlane::Row const& alignPlane,
                         TrkAlignPanel::Row const& alignPanel,
                         Tracker const& nominalTracker, 
                         StrawResponse const& strawRes) {

  std::vector<double> result_locals;
  std::vector<double> result_globals;

  // copy track into this scope - we will modify the params in place a lot
  CosmicTimeTrack track(_track);

  // same here also
  std::vector<double> globals{alignPlane.dx(), alignPlane.dy(), alignPlane.dz(),
                              alignPlane.rx(), alignPlane.ry(), alignPlane.rz(),

                              alignPanel.dx(), alignPanel.dy(), alignPanel.dz(),
                              alignPanel.rx(), alignPanel.ry(), alignPanel.rz()};

  // locals first
  for (size_t paramIdx = 0; paramIdx < track.npars(); ++paramIdx) {
    result_locals.emplace_back(_numericalDerivative(straw, track, globals, nominalTracker, strawRes,
                                                   false, paramIdx, 1e-7));
  }

  for (size_t paramIdx = 0; paramIdx < globals.size(); ++paramIdx) {
    result_globals.emplace_back(
        _numericalDerivative(straw, track, globals, nominalTracker, strawRes, true, paramIdx, 1e-7));
  }

  std::pair<std::vector<double>, std::vector<double>> result{result_locals, result_globals};

  return result;
}

} // namespace AlignmentUtilities
