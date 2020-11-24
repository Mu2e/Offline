// Ryunosuke O'Neil, 2020
// roneil@fnal.gov
// ryunoneil@gmail.com

#include "TrackerAlignment/inc/AlignmentUtilities.hh"
#include "TrackerAlignment/inc/AlignmentDerivatives.hh"
#include "TrackerAlignment/inc/MilleDataWriter.hh"
#include "DataProducts/inc/StrawId.hh"
#include "RecoDataProducts/inc/ComboHit.hh"
#include "TrackerConditions/inc/StrawResponse.hh"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "CLHEP/Vector/ThreeVector.h"
#include "Minuit2/MnUserCovariance.h"

#include <cstdint>
#include <iostream>
#include <iterator>
#include <vector>
using namespace mu2e;

namespace AlignmentUtilities {
  using xyzVec = CLHEP::Hep3Vector; // switch to XYZVec TODO

bool testDerivatives(
    TwoLinePCA const& expected_pca,
    Tracker const& alignedTracker,
    CosmicTimeTrack const& track,
    StrawId const& strawId,
    TrkAlignParams const&rowpl,
    TrkAlignParams const&rowpa,
    Tracker const& nominalTracker,
    StrawResponse const& strawRes) {
  double tolerance = 0.5;

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

  std::vector<double> numLocal, numGlobal, anaLocal, anaGlobal;
  std::tie(anaLocal, anaGlobal) = analyticalDerivatives(track, strawId, rowpl, rowpa, nominalTracker, driftvel);
  std::tie(numLocal, numGlobal) = numericalDerivatives(track, strawId, rowpl, rowpa, nominalTracker, strawRes);

  // compare local derivatives

  if (numLocal.size() != anaLocal.size()) {
    std::cerr << "size mismatch  (local)" << std::endl;
    return false;
  }
  
  for (size_t i = 0; i < numLocal.size(); ++i) {
     if (std::abs(anaLocal[i] - numLocal[i]) > tolerance) {
      std::cerr << "local derivative mismatch(idx " << i << "): diff = " 
                << std::abs(anaLocal[i] - numLocal[i]) 
                << std::endl;
       return false;
     }
  }

  // compare global derivatives
  if (numGlobal.size() != anaGlobal.size()) {
    std::cerr << "size mismatch  (global)" << std::endl;
    return false;
  }
  for (size_t i = 0; i < numGlobal.size(); ++i) {
     if (std::abs(anaGlobal[i] - numGlobal[i]) > tolerance) {
       std::cerr << "global derivative mismatch: diff = " 
                << std::abs(anaGlobal[i] - numGlobal[i]) 
                << std::endl;
       return false;
     }
  }

  return true;
}

std::pair<std::vector<double>, std::vector<double>> 
  analyticalDerivatives(CosmicTimeTrack const& track,
    StrawId const& strawId,
    TrkAlignParams const&rowpl,
    TrkAlignParams const&rowpa,
    Tracker const& nominalTracker,
    double const& driftvel) {

  auto const& plane_origin = nominalTracker.getPlane(strawId).origin();
  auto const& panel_origin = nominalTracker.getPanel(strawId).origin();
  auto const& nominalStraw = nominalTracker.getStraw(strawId);
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
std::pair<Hep3Vector, Hep3Vector> alignStraw(Tracker const& tracker, 
                                             StrawId const& strawId,
                                             TrkAlignParams const& align_tracker,
                                             TrkAlignParams const& align_plane,
                                             TrkAlignParams const& align_panel,
					     TrkStrawEndAlign const& align_straw) {  

  auto const& plane = tracker.getPlane(strawId);
  auto const& panel = tracker.getPanel(strawId);
  auto const& straw = tracker.getStraw(strawId);

  // nominal place of plane in the tracker
  auto plane_to_ds = plane.planeToDS();
  // nominal panel in the plane
  auto panel_to_plane = plane.dsToPlane()*panel.panelToDS();
  // chained alignment including panel
  auto aligned_panel_to_ds = align_tracker.transform() * (plane_to_ds * align_plane.transform()) * (panel_to_plane * align_panel.transform());
  // Move wire back to the nominal panel frame (UVW)
  auto ds_to_panel = panel.dsToPanel();
  std::array<xyzVec,2> wireends;
  for(int iend=0;iend < StrawEnd::nends; iend++){
    auto end = static_cast<StrawEnd::End>(iend);
    // include the straw end alignment
    auto wireend_UVW = ds_to_panel*straw.wireEnd(end) + align_straw.wireDeltaUVW(end);
    wireends[iend] = aligned_panel_to_ds*wireend_UVW;
  }
// compute the net position and direction from these
  auto aligned_strawmid = 0.5*(wireends[0] + wireends[1]); 
  auto aligned_strawdir = (wireends[StrawEnd::cal]- wireends[StrawEnd::hv]).unit(); // convention  is direction points from HV to Cal
  return {aligned_strawmid, aligned_strawdir};
}



/* Numerical partial DOCA/TOCA derivatives
 *
 */ 


namespace {

// do the alignment parameters really need to be double?
  double tocaGlobalDep(CosmicTimeTrack const& track, StrawId const& strawId,
                                std::vector<double> const& globals, Tracker const& nominalTracker, StrawResponse const& strawRes) {


    TrkAlignParams align_tracker{strawId, StrawIdMask::tracker, 0, 0, 0, 0, 0, 0};
    TrkAlignParams align_plane{strawId, StrawIdMask::plane, globals[0], globals[1], globals[2], globals[3], globals[4], globals[5]};
    TrkAlignParams align_panel{strawId, StrawIdMask::panel, globals[6], globals[7],  globals[8], globals[9], globals[10], globals[11]};
    TrkStrawEndAlign align_straw{strawId.uniqueStraw(), strawId, 0,0,0,0,0,0,0,0}; // not sure how to really initialize this FIXME!

    // returns pair of vectors { straw_pos, straw_dir }
    Hep3Vector straw_pos, straw_dir;
    std::tie(straw_pos, straw_dir) = alignStraw(nominalTracker, 
                                                strawId, align_tracker, align_plane, align_panel, align_straw);

    TwoLinePCA pca(track.intercept(), track.direction(), straw_pos, straw_dir);

    double traj_time = (pca.point1() - track.intercept()).dot(track.direction()) / 299.9;
    double d2t_doca = strawRes.driftDistanceToTime(strawId, pca.dca(), 0);
    double t_offset = strawRes.driftTimeOffset(strawId, 0, 0, pca.dca());

    double predictedTime = traj_time + t_offset + track.params[CosmicTimeTrack::t0] + d2t_doca;

    return predictedTime;
  }

  double docaGlobalDep(CosmicTimeTrack const& track, StrawId const& strawId,
                                std::vector<double> const& globals, Tracker const& nominalTracker, StrawResponse const& strawRes) {
// this function should be consolidated with the above FIXME!
    TrkAlignParams align_tracker{strawId, StrawIdMask::tracker, 0, 0, 0, 0, 0, 0};
    TrkAlignParams align_plane{strawId, StrawIdMask::plane, globals[0], globals[1], globals[2], globals[3], globals[4], globals[5]};
    TrkAlignParams align_panel{strawId, StrawIdMask::panel, globals[6], globals[7],  globals[8], globals[9], globals[10], globals[11]}; 
    TrkStrawEndAlign align_straw{strawId.uniqueStraw(), strawId, 0,0,0,0,0,0,0,0}; // not sure how to really initialize this FIXME!
    // returns pair of vectors { straw_pos, straw_dir }
    Hep3Vector straw_pos, straw_dir;
    std::tie(straw_pos, straw_dir) = alignStraw(nominalTracker, 
                                                strawId, align_tracker, align_plane, align_panel,align_straw);

    TwoLinePCA pca(track.intercept(), track.direction(), straw_pos, straw_dir);

    int ambig = hitAmbiguity(track, straw_pos, straw_dir);
    double doca = ambig * pca.dca();

    return doca;
  }


  // not meant to be called from outside of this namespace
  double _numericalDerivative(StrawId const& straw, CosmicTimeTrack& track,
                            std::vector<double>& globals, Tracker const& nominalTracker,
                            StrawResponse const& strawRes, bool isGlobalParam,
                            size_t const& paramIdx, double step_size, bool useTimeDomain) {
    // calculate numerical partial derivative wrt param at paramIdx in either
    // local, or global param array

    double x;

    if (isGlobalParam) {
      x = globals[paramIdx];
      globals[paramIdx] = x + step_size;
    } else {
      x = track.params[paramIdx];
      track.params[paramIdx] = x + step_size;
    }

    double pdiff;
    
    if (useTimeDomain) {
      pdiff = tocaGlobalDep(track, straw, globals, nominalTracker, strawRes);
    }
    else {
      pdiff = docaGlobalDep(track, straw, globals, nominalTracker, strawRes);
    }

    if (isGlobalParam) {
      globals[paramIdx] = x - step_size;
    } else {
      track.params[paramIdx] = x - step_size;
    }

    if (useTimeDomain) {
      pdiff -= tocaGlobalDep(track, straw, globals, nominalTracker, strawRes);
    }
    else {
      pdiff -= docaGlobalDep(track, straw, globals, nominalTracker, strawRes);
    }

    pdiff /= (2.0 * step_size);

    if (isGlobalParam) {
      globals[paramIdx] = x;
    } else {
      track.params[paramIdx] = x;
    }

    return pdiff;
  }
}

std::pair<std::vector<double>, std::vector<double>>
numericalDerivatives(CosmicTimeTrack const& _track, StrawId const& straw,
                         TrkAlignParams const& alignPlane,
                         TrkAlignParams const& alignPanel,
                         Tracker const& nominalTracker, 
                         StrawResponse const& strawRes,
                         bool useTimeDomain) {

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
    result_locals.emplace_back(
        _numericalDerivative(straw, track, globals, nominalTracker, strawRes,
                              false, paramIdx, 1e-7, useTimeDomain));
  }

  for (size_t paramIdx = 0; paramIdx < globals.size(); ++paramIdx) {
    result_globals.emplace_back(
        _numericalDerivative(straw, track, globals, nominalTracker, strawRes, 
                              true, paramIdx, 1e-7, useTimeDomain));
  }


  return {result_locals, result_globals};
}

} // namespace AlignmentUtilities
