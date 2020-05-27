// Ryunosuke O'Neil, 2019
// roneil@fnal.gov
// ryunoneil@gmail.com

// A module to collect Cosmic NoField tracks and write out 'Mille' data files used as input to a
// Millepede-II alignment fit.

// Consult README.md for more information

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <fstream>
#include <iostream>
#include <memory>
#include <ostream>
#include <stdint.h>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "CosmicReco/inc/PDFFit.hh"
#include "DbTables/inc/TrkAlignPanel.hh"
#include "DbTables/inc/TrkAlignPlane.hh"
#include "GeneralUtilities/inc/BitMap.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Minuit2/MnUserCovariance.h"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

#include "boost/math/distributions/chi_squared.hpp"
#include "boost/math/distributions/normal.hpp"

#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Core/ProducerTable.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"
#include "art_root_io/TFileService.h"

#include "TrackerConditions/inc/StrawResponse.hh"
#include "TrackerGeom/inc/Panel.hh"
#include "TrackerGeom/inc/Plane.hh"
#include "TrackerGeom/inc/Straw.hh"
#include "TrackerGeom/inc/Tracker.hh"

#include "DbService/inc/DbHandle.hh"
#include "ProditionsService/inc/ProditionsHandle.hh"

#include "RecoDataProducts/inc/ComboHit.hh"
#include "RecoDataProducts/inc/CosmicTrack.hh"
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"

#include "DataProducts/inc/StrawId.hh"
#include "DataProducts/inc/XYZVec.hh"

#include "TAxis.h"
#include "TH1F.h"
#include "TMatrixDSym.h"
#include "TTree.h"

#include "canvas/Utilities/Exception.h"
#include "canvas/Utilities/InputTag.h"

#include "CLHEP/Vector/ThreeVector.h"

#include "cetlib_except/exception.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Comment.h"
#include "fhiclcpp/types/Name.h"
#include "fhiclcpp/types/Table.h"
#include "fhiclcpp/types/detail/validationException.h"

#include "TrackerAlignment/inc/AlignmentDerivatives.hh"
#include "TrackerAlignment/inc/Mille.h"

namespace art {
class Run;
} // namespace art

using namespace mu2e;
namespace mu2e {

class AlignTrackCollector : public art::EDAnalyzer {
private:
  static constexpr int MAX_NHITS = 100;

  // Tree and tree fill members
  TTree* diagtree;

  Int_t nHits;
  Float_t doca_residual[MAX_NHITS];
  Float_t time_residual[MAX_NHITS];
  Float_t residual_err[MAX_NHITS];
  Float_t doca_resid_err[MAX_NHITS];
  Float_t drift_reso[MAX_NHITS];

  Float_t pull_doca[MAX_NHITS];
  Float_t pull_hittime[MAX_NHITS];

  Float_t doca[MAX_NHITS];
  Float_t time[MAX_NHITS];
  Int_t plane_uid[MAX_NHITS];
  Int_t panel_uid[MAX_NHITS];

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

public:
  size_t _dof_per_plane = 6; // dx, dy, dz, a, b, g (translation, rotation)
  size_t _dof_per_panel = 6; // dx, dy, dz, a, b, g (translation, rotation)
  size_t _ndof = StrawId::_nplanes * _dof_per_plane + StrawId::_nupanels * _dof_per_panel;
  size_t _expected_dofs = _dof_per_panel + _dof_per_plane;

  struct Config {
    using Name = fhicl::Name;
    using Comment = fhicl::Comment;

    fhicl::Atom<int> diaglvl{Name("diagLevel"), Comment("diagnostic level")};

    fhicl::Atom<art::InputTag> costag{Name("CosmicTrackSeedCollection"),
                                      Comment("tag for cosmic track seed collection")};

    fhicl::Atom<std::string> millefile{Name("TrackDataOutputFile"),
                                       Comment("Output filename for Millepede track data file")};

    fhicl::Atom<std::string> labelsfile{
        Name("LabelsOutputFile"), Comment("Output filename for Millepede label ID's (debug only)"),
        "none"};

    fhicl::Atom<std::string> tracktype{
        Name("TrackType"),
        Comment("The type of track to collect. Default: CosmicTrackSeedCollection"),
        "CosmicTrackSeedCollection"};

    fhicl::Atom<int> minplanetraverse{Name("MinTraversedPlanes"),
                                      Comment("How many planes must be traversed for a track "
                                              "to be accepted. 0: does not apply the cut."),
                                      3};

    fhicl::Atom<int> minpaneltraverse{
        Name("MinTraversedPanelsPerPlane"),
        Comment("How many panels must be traversed PER PLANE. 0: does not apply the cut."), 0};

    fhicl::Atom<double> maxpvalue{Name("MaxPValue"),
                                  Comment("Require that the track p-value < MaxPValue"), 1};

    fhicl::Atom<double> maxtimeres{Name("MaxTimeRes"),
                                   Comment("Require that the maximum ABSOLUTE time residual over "
                                           "all track hits < MaxTimeRes. Setting a "
                                           "negative value does not apply the cut."),
                                   -1.0};

    fhicl::Atom<int> mintrackhits{
        Name("MinTrackSH"), Comment("Require that the minimum straw hits in a track > MinTrackSH."),
        0};

    fhicl::Atom<bool> usetimeresid{
        Name("UseTimeDomain"),
        Comment("Write the alignment data in the time domain. i.e. measurement = Time "
                "Residual, sigma = driftTimeError."),
        true};

    fhicl::Atom<bool> nopaneldofs{Name("NoPanelDOFs"), Comment("remove panel DOFs"), true};

    fhicl::Atom<bool> noplanerotations{Name("NoPlaneRotations"),
                                       Comment("Remove Plane rotation DOFs"), true};

    fhicl::Atom<bool> useplanefilter{
        Name("PlaneFilter"),
        Comment("Only write hit measurements made in specified planes (PlaneFilterList)"), false};

    fhicl::Sequence<int> planefilterlist{
        Name("PlaneFilterList"),
        Comment("Only write measurements made in these planes (int list [0,35])"),
    };
  };

  typedef art::EDAnalyzer::Table<Config> Parameters;

  void beginJob();
  void endJob();
  void beginRun(art::Run const&);
  void analyze(art::Event const&);
  bool filter_CosmicTrackSeedCollection(art::Event const& event, Tracker const& tracker,
                                        Tracker const& nominalTracker, StrawResponse const& _srep,
                                        CosmicTrackSeedCollection const& _coscol);

  int getLabel(int const&, int const&, int const&);
  std::vector<int> generateDOFLabels(uint16_t plane, uint16_t panel);
  std::vector<int> generateDOFLabels(StrawId const& strw);
  void writeMillepedeParams(TrkAlignPlane const& alignConstPlanes,
                            TrkAlignPanel const& alignConstPanels);

  AlignTrackCollector(const Parameters& conf) :
      art::EDAnalyzer(conf), _diag(conf().diaglvl()), _costag(conf().costag()),
      _output_filename(conf().millefile()), _labels_filename(conf().labelsfile()),
      track_type(conf().tracktype()), min_plane_traverse(conf().minplanetraverse()),
      min_panel_traverse_per_plane(conf().minpaneltraverse()), max_pvalue(conf().maxpvalue()),
      max_timeres(conf().maxtimeres()), min_track_hits(conf().mintrackhits()),
      use_timeresid(conf().usetimeresid()), no_panel_dofs(conf().nopaneldofs()),
      no_plane_rotations(conf().noplanerotations()), use_plane_filter(conf().useplanefilter()),
      plane_filter_list(conf().planefilterlist()),

      wroteMillepedeParams(false) {

    if (no_panel_dofs) {
      _dof_per_panel = 0;
    }

    if (no_plane_rotations) {
      _dof_per_plane = 3;
    }

    _ndof = StrawId::_nplanes * _dof_per_plane + StrawId::_nupanels * _dof_per_panel;
    _expected_dofs = _dof_per_panel + _dof_per_plane;

    if (_diag > 0) {
      std::cout << "AlignTrackCollector: Total number of plane degrees of freedom = "
                << StrawId::_nplanes * _dof_per_plane << std::endl;

      std::cout << "AlignTrackCollector: Total number of panel degrees of freedom = "
                << StrawId::_nupanels * _dof_per_panel << std::endl;
    }
  }

  virtual ~AlignTrackCollector() {}

  Config _conf;

  int _diag;
  art::InputTag _costag;
  std::string _output_filename;
  std::string _labels_filename;
  std::string track_type;

  int min_plane_traverse;
  int min_panel_traverse_per_plane;
  double max_pvalue;
  double max_timeres;
  int min_track_hits;
  bool use_timeresid;
  bool no_panel_dofs;
  bool no_plane_rotations;

  bool use_plane_filter;
  std::vector<int> plane_filter_list;

  bool use_db;

  bool wroteMillepedeParams;

  std::unique_ptr<Mille> millepede;
  const CosmicTrackSeedCollection* _coscol;
  const Tracker* _tracker;

  size_t tracks_written = 0;

  ProditionsHandle<Tracker> _proditionsTracker_h;
  ProditionsHandle<StrawResponse> srep_h;

  std::unique_ptr<DbHandle<TrkAlignPlane>> _trkAlignPlane_h;
  std::unique_ptr<DbHandle<TrkAlignPanel>> _trkAlignPanel_h;
};

void AlignTrackCollector::beginJob() {
  millepede = std::make_unique<Mille>(_output_filename.c_str());
  _trkAlignPlane_h = std::make_unique<DbHandle<TrkAlignPlane>>();
  _trkAlignPanel_h = std::make_unique<DbHandle<TrkAlignPanel>>();

  if (_diag > 0) {
    art::ServiceHandle<art::TFileService> tfs;

    diagtree = tfs->make<TTree>("tracks", "Tracks collected for an alignment iteration");
    diagtree->Branch("nHits", &nHits, "nHits/I");
    diagtree->Branch("doca_resid", &doca_residual, "doca_resid[nHits]/F");
    diagtree->Branch("time_resid", &time_residual, "time_resid[nHits]/F");
    diagtree->Branch("doca_resid_err", &doca_resid_err, "doca_resid_err[nHits]/F");
    diagtree->Branch("drift_res", &drift_reso, "drift_res[nHits]/F");

    diagtree->Branch("resid_err", &residual_err, "resid_err[nHits]/F");

    diagtree->Branch("pull_doca", &pull_doca, "pull_doca[nHits]/F");
    diagtree->Branch("pull_hittime", &pull_hittime, "pull_doca[nHits]/F");

    diagtree->Branch("doca", &doca, "doca[nHits]/F");
    diagtree->Branch("time", &time, "time[nHits]/F");
    diagtree->Branch("plane", &plane_uid, "plane[nHits]/I");
    diagtree->Branch("panel", &panel_uid, "panel[nHits]/I");

    diagtree->Branch("A0", &A0, "A0/D");
    diagtree->Branch("A1", &A1, "A1/D");
    diagtree->Branch("B0", &B0, "B0/D");
    diagtree->Branch("B1", &B1, "B1/D");
    diagtree->Branch("T0", &T0, "T0/D");

    diagtree->Branch("chisq", &chisq, "chisq/D");
    diagtree->Branch("chisq_doca", &chisq_doca, "chisq_doca/D");

    diagtree->Branch("ndof", &ndof, "ndof/I");
    diagtree->Branch("pvalue", &pvalue, "pvalue/D");

    diagtree->Branch("panels_trav", &panels_trav, "panels_trav/I");
    diagtree->Branch("planes_trav", &planes_trav, "planes_trav/I");
  }
}

void AlignTrackCollector::writeMillepedeParams(TrkAlignPlane const& alignConstPlanes,
                                               TrkAlignPanel const& alignConstPanels) {
  // write a params.txt telling millepede what the alignment constants
  // were for this track collection iteration

  std::ofstream output_file("mp_params.txt"); // FIXME: fcl param

  output_file << "Parameter" << std::endl;

  for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
    auto const& rowpl = alignConstPlanes.rowAt(p);

    if (use_plane_filter && std::find(plane_filter_list.begin(), plane_filter_list.end(), p) ==
                                plane_filter_list.end()) {
      // we're not interested in this plane!
      continue;
    }

    std::vector<float> pl_consts{rowpl.dx(), rowpl.dy(), rowpl.dz(),
                                 rowpl.rx(), rowpl.ry(), rowpl.rz()};

    for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
      // label(int)   initial value (float)   presigma (-ve: fixed, 0: variable)
      output_file << getLabel(1, p, dof_n) << "  " << pl_consts[dof_n] << "  0" << std::endl;
    }
  }

  output_file << std::endl << "! Panel DOFs follow:" << std::endl;

  for (uint16_t p = 0; p < StrawId::_nupanels; ++p) {
    auto const& rowpa = alignConstPanels.rowAt(p);
    int plane = p % 6;

    if (use_plane_filter && std::find(plane_filter_list.begin(), plane_filter_list.end(), plane) ==
                                plane_filter_list.end()) {
      continue;
    }

    std::vector<float> pa_consts{rowpa.dx(), rowpa.dy(), rowpa.dz(),
                                 rowpa.rx(), rowpa.ry(), rowpa.rz()};

    for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
      // label(int)   initial value (float)   presigma (-ve: fixed, 0: variable)
      output_file << getLabel(2, p, dof_n) << "  " << pa_consts[dof_n] << "  0" << std::endl;
    }
  }

  output_file.close();
}

std::vector<int> AlignTrackCollector::generateDOFLabels(StrawId const& strw) {
  return generateDOFLabels(strw.getPlane(), strw.getPanel());
}

std::vector<int> AlignTrackCollector::generateDOFLabels(uint16_t plane, uint16_t panel) {
  std::vector<int> labels;
  labels.reserve(_dof_per_plane + _dof_per_panel);

  for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
    labels.push_back(getLabel(1, plane, dof_n));
  }

  if (!no_panel_dofs) {
    for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
      labels.push_back(getLabel(2, panel, dof_n));
    }
  }
  return labels; // should 'move' the vector, not copy it
}

void AlignTrackCollector::beginRun(art::Run const& run) {
  GeomHandle<Tracker> track;
  _tracker = track.get();

  return;
}

void AlignTrackCollector::endJob() {
  // ensure the binary Mille file is explicitly closed in endJob,
  // not just when this module goes out of scope.
  // File close is handled in Mille destructor
  millepede->~Mille();

  if (_diag > 0) {
    std::cout << "AlignTrackCollector: wrote " << tracks_written << " tracks to "
              << _output_filename << std::endl;
  }
}

bool AlignTrackCollector::filter_CosmicTrackSeedCollection(
    art::Event const& event, Tracker const& tracker, Tracker const& nominalTracker,
    StrawResponse const& _srep, CosmicTrackSeedCollection const& coscol) {

  // get alignment parameters for this event
  // N.B. alignment parameters MUST be unchanged for the entire job...
  // FIXME! check and enforce above
  auto alignConsts_planes = _trkAlignPlane_h->get(event.id());
  auto alignConsts_panels = _trkAlignPanel_h->get(event.id());

  if (!wroteMillepedeParams) {
    writeMillepedeParams(alignConsts_planes, alignConsts_panels);
    wroteMillepedeParams = true;
  }

  bool wrote_track = false; // did we write any tracks at all?

  // dedicated to CosmicTrackSeedCollection
  for (CosmicTrackSeed const& sts : coscol) {
    CosmicTrack const& st = sts._track;
    TrkFitFlag const& status = sts._status;

    if (!status.hasAllProperties(TrkFitFlag::helixOK)) {
      continue;
    }

    if (!st.converged || !st.minuit_converged) {
      continue;
    }

    if (isnan(st.MinuitParams.A0)) {
      continue;
    }

    std::set<uint16_t> planes_traversed;
    std::set<uint16_t> panels_traversed;

    A0 = st.MinuitParams.A0;
    A1 = st.MinuitParams.A1;
    B0 = st.MinuitParams.B0;
    B1 = st.MinuitParams.B1;
    T0 = st.MinuitParams.T0;

    CLHEP::Hep3Vector intercept(A0, 0, B0);
    CLHEP::Hep3Vector dir(A1, -1, B1);
    dir = dir.unit();

    GaussianDriftFit fit_object(sts._straw_chits, _srep, &tracker);

    chisq = 0;
    chisq_doca = 0;
    ndof = 0;
    pvalue = 0;
    nHits = 0;

    // for the max timeresidual track quality cut
    double max_time_res_track = -1;

    bool wrote_hits = false; // did we write any hits for the track?
    bool bad_track = false;

    // V = cov( measurements )
    TMatrixD meas_cov(sts._straw_chits.size(), sts._straw_chits.size());
    meas_cov.Zero();

    // H = partial d (residuals)/ d (track params)
    TMatrixD resid_local_derivs(sts._straw_chits.size(), 5);
    resid_local_derivs.Zero();

    // C = cov ( track parameters )
    // FIXME! Minuit has its own conventions for storing
    // covariance matrices..!
    ROOT::Minuit2::MnUserCovariance cov(sts._track.MinuitParams.cov, 5);
    TMatrixD track_cov(5, 5);
    for (size_t r = 0; r < 5; r++) {
      for (size_t c = 0; c < 5; c++) {
        track_cov(r, c) = cov(r, c);
      }
    }

    // TODO: Check positive semi-definiteness

    std::cout << "Constructed track covariance...?" << std::endl;

    std::vector<float> residuals;
    std::vector<std::vector<float>> global_derivs_temp;
    std::vector<std::vector<float>> local_derivs_temp;
    std::vector<std::vector<int>> labels_temp;

    // get residuals and their derivatives with respect
    // to all local and global parameters
    // get also plane id hit by straw hits
    for (ComboHit const& straw_hit : sts._straw_chits) {
      // straw and plane info
      StrawId const& straw_id = straw_hit.strawId();
      Straw const& straw = tracker.getStraw(straw_id);
      Straw const& nominalStraw = nominalTracker.getStraw(straw_id);

      auto plane_id = straw_id.getPlane();
      auto panel_uuid = straw_id.uniquePanel();
      auto panel_id = straw_id.getPanelId();

      // geometry info for derivatives
      auto const& plane_origin = nominalTracker.getPlane(plane_id).origin();
      auto const& panel_origin = nominalTracker.getPanel(panel_id).straw0MidPoint();
      auto const& straw_mp = nominalStraw.getMidPoint();
      auto const& wire_dir = nominalStraw.getDirection();
      auto const& rowpl = alignConsts_planes.rowAt(plane_id);
      auto const& rowpa = alignConsts_panels.rowAt(panel_uuid);

      TwoLinePCA pca(straw.getMidPoint(), straw.getDirection(), intercept, dir);

      double driftvel = _srep.driftInstantSpeed(straw_id, pca.dca(), 0);

      // now calculate the derivatives
      auto derivativesLocal = CosmicTrack_DCA_LocalDeriv(
          A0, B0, A1, B1, T0, rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(), rowpl.ry(),
          rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(), rowpa.ry(), rowpa.rz(),

          straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(), wire_dir.y(), wire_dir.z(),
          plane_origin.x(), plane_origin.y(), plane_origin.z(), panel_origin.x(), panel_origin.y(),
          panel_origin.z(), driftvel);

      auto derivativesGlobal = CosmicTrack_DCA_GlobalDeriv(
          A0, B0, A1, B1, T0, rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(), rowpl.ry(),
          rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(), rowpa.ry(), rowpa.rz(),

          straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(), wire_dir.y(), wire_dir.z(),

          plane_origin.x(), plane_origin.y(), plane_origin.z(), panel_origin.x(), panel_origin.y(),
          panel_origin.z(), driftvel);

      double dca_resid = fit_object.DOCAresidual(straw_hit, sts);
      double drift_res_dca = _srep.driftDistanceError(straw_hit.strawId(), 0, 0, pca.dca());
      double signdca = (pca.s2() > 0 ? pca.dca() : -pca.dca());

      double time_resid = fit_object.TimeResidual(straw_hit, sts);
      double drift_res = _srep.driftTimeError(straw_hit.strawId(), 0, 0, pca.dca());

      // The following are based on reco performed using current alignment parameters
      // FIXME: confusing



      // FIXME! use newly implemented chisq function in fit object
      chisq += pow(time_resid / drift_res, 2);
      chisq_doca += pow(dca_resid / drift_res_dca, 2);

      if (isnan(dca_resid) || isnan(time_resid) || isnan(drift_res)) {
        bad_track = true;
        continue;
      }

      planes_traversed.insert(plane_id);
      panels_traversed.insert(panel_uuid);

      doca_residual[nHits] = dca_resid;
      time_residual[nHits] = time_resid;
      doca_resid_err[nHits] = drift_res_dca;
      pull_doca[nHits] = dca_resid / drift_res_dca;
      pull_hittime[nHits] = time_resid / drift_res;
      drift_reso[nHits] = drift_res;
      doca[nHits] = pca.dca();
      time[nHits] = straw_hit.time();
      panel_uid[nHits] = panel_uuid;
      plane_uid[nHits] = plane_id;

      if (_diag > 1) {
        std::cout << "pl" << plane_id << " pa" << panel_uuid << ": dcaresid " << dca_resid << " +- "
                  << drift_res_dca << std::endl;
        std::cout << ": timeresid " << time_resid << " +- " << drift_res << std::endl;
      }

      // diagnostics!
      // diagnostics!
      // diagnostics!
      // diagnostics!
      // diagnostics!
      // diagnostics!

      if (_diag > 4) {
        // FIXME!
        // move to another place
        double generated_doca = CosmicTrack_DCA(
            A0, B0, A1, B1, T0, rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(), rowpl.ry(),
            rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(), rowpa.ry(), rowpa.rz(),

            straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(), wire_dir.y(), wire_dir.z(),
            plane_origin.x(), plane_origin.y(), plane_origin.z(), panel_origin.x(),
            panel_origin.y(), panel_origin.z(), driftvel);

        double diff = std::abs(signdca - generated_doca);
        std::cout << "doca: " << signdca << ", gendoca: " << generated_doca << ", diff: " << diff
                  << ", s1: " << pca.s1() << ", s2: " << pca.s2() << std::endl;

        if (diff > 1e-10) {
          std::cout << "----------------------------------" << std::endl;
          std::cout
              << "WARNING! Output of generated DOCA function (after alignment) not consistent!"
              << std::endl;
          std::cout << "----------------------------------" << std::endl;

          // throw cet::exception("ALIGNMENT") << "Output of generated functions
          // (AlignmentDerivatives) are not consistent!";
        }

        std::cout << "dr/d([A0,B0,A1,B1,T0]): [" << derivativesLocal[0] // A0
                  << ", " << derivativesLocal[1]                        // B0
                  << ", " << derivativesLocal[2]                        // A1
                  << ", " << derivativesLocal[3]                        // B1
                  << ", " << derivativesLocal[4];                       // T0
        std::cout << "]" << std::endl;

        std::cout << "dr/d(plane" << plane_id << "[x, y, z]): [" << derivativesGlobal[0] // x
                  << ", " << derivativesGlobal[1]                                        // y
                  << ", " << derivativesGlobal[2];                                       // z
        std::cout << "]" << std::endl;

        // quick and dirty numerical derivative estimation
        // This is only performed using CosmicTrack_DCA once confirmed to be consistent with
        // TwoLinePCA
        // TODO: move to a utility class?
        // TODO: avoid DRY problems

        derivativesLocal.clear();
        derivativesGlobal.clear();

        double h = 1e-7;

        // double linvel = 0.0625;

        // PARTIAL DOCA DERIVATIVE: A0

        double diff_a = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0 + h, B0, A1, B1, T0, rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(),
                            rowpl.ry(), rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(),
                            rowpa.ry(), rowpa.rz(),

                            straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(), wire_dir.y(),
                            wire_dir.z(), plane_origin.x(), plane_origin.y(), plane_origin.z(),
                            panel_origin.x(), panel_origin.y(), panel_origin.z(), driftvel),
            0);

        double diff_b = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0 - h, B0, A1, B1, T0, rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(),
                            rowpl.ry(), rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(),
                            rowpa.ry(), rowpa.rz(),

                            straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(), wire_dir.y(),
                            wire_dir.z(), plane_origin.x(), plane_origin.y(), plane_origin.z(),
                            panel_origin.x(), panel_origin.y(), panel_origin.z(), driftvel),
            0);

        diff = (diff_a - diff_b) / (2.0 * h);
        std::cout << "numerical dr/d(A0) = " << diff << std::endl;

        derivativesLocal.push_back(diff);

        // PARTIAL DOCA DERIVATIVE: B0

        diff_a = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0, B0 + h, A1, B1, T0, rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(),
                            rowpl.ry(), rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(),
                            rowpa.ry(), rowpa.rz(),

                            straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(), wire_dir.y(),
                            wire_dir.z(), plane_origin.x(), plane_origin.y(), plane_origin.z(),
                            panel_origin.x(), panel_origin.y(), panel_origin.z(), driftvel),
            0);

        diff_b = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0, B0 - h, A1, B1, T0, rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(),
                            rowpl.ry(), rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(),
                            rowpa.ry(), rowpa.rz(),

                            straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(), wire_dir.y(),
                            wire_dir.z(), plane_origin.x(), plane_origin.y(), plane_origin.z(),
                            panel_origin.x(), panel_origin.y(), panel_origin.z(), driftvel),
            0);

        diff = (diff_a - diff_b) / (2.0 * h);
        std::cout << "numerical dr/d(B0) = " << diff << std::endl;
        derivativesLocal.push_back(diff);

        // PARTIAL DOCA DERIVATIVE: A1

        diff_a = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0, B0, A1 + h, B1, T0, rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(),
                            rowpl.ry(), rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(),
                            rowpa.ry(), rowpa.rz(),

                            straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(), wire_dir.y(),
                            wire_dir.z(), plane_origin.x(), plane_origin.y(), plane_origin.z(),
                            panel_origin.x(), panel_origin.y(), panel_origin.z(), driftvel),
            0);

        diff_b = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0, B0, A1 - h, B1, T0, rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(),
                            rowpl.ry(), rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(),
                            rowpa.ry(), rowpa.rz(),

                            straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(), wire_dir.y(),
                            wire_dir.z(), plane_origin.x(), plane_origin.y(), plane_origin.z(),
                            panel_origin.x(), panel_origin.y(), panel_origin.z(), driftvel),
            0);

        diff = (diff_a - diff_b) / (2.0 * h);
        std::cout << "numerical dr/d(A1) = " << diff << std::endl;
        derivativesLocal.push_back(diff);

        // PARTIAL DOCA DERIVATIVE: B1

        diff_a = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0, B0, A1, B1 + h, T0, rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(),
                            rowpl.ry(), rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(),
                            rowpa.ry(), rowpa.rz(), straw_mp.x(), straw_mp.y(), straw_mp.z(),
                            wire_dir.x(), wire_dir.y(), wire_dir.z(), plane_origin.x(),
                            plane_origin.y(), plane_origin.z(), panel_origin.x(), panel_origin.y(),
                            panel_origin.z(), driftvel),
            0);

        diff_b = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0, B0, A1, B1 - h, T0, rowpl.dx(), rowpl.dy(), rowpl.dz(), rowpl.rx(),
                            rowpl.ry(), rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(),
                            rowpa.ry(), rowpa.rz(), straw_mp.x(), straw_mp.y(), straw_mp.z(),
                            wire_dir.x(), wire_dir.y(), wire_dir.z(), plane_origin.x(),
                            plane_origin.y(), plane_origin.z(), panel_origin.x(), panel_origin.y(),
                            panel_origin.z(), driftvel),
            0);

        diff = (diff_a - diff_b) / (2.0 * h);
        std::cout << "numerical dr/d(B1) = " << diff << std::endl;
        derivativesLocal.push_back(diff);

        // T0
        derivativesLocal.push_back(-1);

        // PARTIAL DOCA DERIVATIVE: PLANE X SHIFT

        diff_a = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0, B0, A1, B1, T0,

                            rowpl.dx() + h,

                            rowpl.dy(), rowpl.dz(), rowpl.rx(), rowpl.ry(), rowpl.rz(), rowpa.dx(),
                            rowpa.dy(), rowpa.dz(), rowpa.rx(), rowpa.ry(), rowpa.rz(),
                            straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(), wire_dir.y(),
                            wire_dir.z(), plane_origin.x(), plane_origin.y(), plane_origin.z(),
                            panel_origin.x(), panel_origin.y(), panel_origin.z(), driftvel),
            0);

        diff_b = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0, B0, A1, B1, T0,

                            rowpl.dx() - h,

                            rowpl.dy(), rowpl.dz(), rowpl.rx(), rowpl.ry(), rowpl.rz(), rowpa.dx(),
                            rowpa.dy(), rowpa.dz(), rowpa.rx(), rowpa.ry(), rowpa.rz(),
                            straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(), wire_dir.y(),
                            wire_dir.z(), plane_origin.x(), plane_origin.y(), plane_origin.z(),
                            panel_origin.x(), panel_origin.y(), panel_origin.z(), driftvel),
            0);

        diff = (diff_a - diff_b) / (2.0 * h);
        std::cout << "numerical dr/d([plane " << plane_id << "]DX) = " << diff << std::endl;
        derivativesGlobal.push_back(diff);

        // PARTIAL DOCA DERIVATIVE: PLANE Y SHIFT

        diff_a = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0, B0, A1, B1, T0, rowpl.dx(), rowpl.dy() + h, rowpl.dz(), rowpl.rx(),
                            rowpl.ry(), rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(),
                            rowpa.ry(), rowpa.rz(), straw_mp.x(), straw_mp.y(), straw_mp.z(),
                            wire_dir.x(), wire_dir.y(), wire_dir.z(), plane_origin.x(),
                            plane_origin.y(), plane_origin.z(), panel_origin.x(), panel_origin.y(),
                            panel_origin.z(), driftvel),
            0);

        diff_b = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0, B0, A1, B1, T0, rowpl.dx(), rowpl.dy() - h, rowpl.dz(), rowpl.rx(),
                            rowpl.ry(), rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(),
                            rowpa.ry(), rowpa.rz(), straw_mp.x(), straw_mp.y(), straw_mp.z(),
                            wire_dir.x(), wire_dir.y(), wire_dir.z(), plane_origin.x(),
                            plane_origin.y(), plane_origin.z(), panel_origin.x(), panel_origin.y(),
                            panel_origin.z(), driftvel),
            0);

        diff = (diff_a - diff_b) / (2.0 * h);
        std::cout << "numerical dr/d([plane " << plane_id << "]DY) = " << diff << std::endl;
        derivativesGlobal.push_back(diff);

        // PARTIAL DOCA DERIVATIVE: PLANE Z SHIFT

        diff_a = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0, B0, A1, B1, T0, rowpl.dx(), rowpl.dy(), rowpl.dz() + h, rowpl.rx(),
                            rowpl.ry(), rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(),
                            rowpa.ry(), rowpa.rz(), straw_mp.x(), straw_mp.y(), straw_mp.z(),
                            wire_dir.x(), wire_dir.y(), wire_dir.z(), plane_origin.x(),
                            plane_origin.y(), plane_origin.z(), panel_origin.x(), panel_origin.y(),
                            panel_origin.z(), driftvel),
            0);

        diff_b = _srep.driftDistanceToTime(
            straw_id,
            CosmicTrack_DCA(A0, B0, A1, B1, T0, rowpl.dx(), rowpl.dy(), rowpl.dz() - h, rowpl.rx(),
                            rowpl.ry(), rowpl.rz(), rowpa.dx(), rowpa.dy(), rowpa.dz(), rowpa.rx(),
                            rowpa.ry(), rowpa.rz(), straw_mp.x(), straw_mp.y(), straw_mp.z(),
                            wire_dir.x(), wire_dir.y(), wire_dir.z(), plane_origin.x(),
                            plane_origin.y(), plane_origin.z(), panel_origin.x(), panel_origin.y(),
                            panel_origin.z(), driftvel),
            0);

        diff = (diff_a - diff_b) / (2.0 * h);
        std::cout << "numerical dr/d([plane " << plane_id << "]DZ) = " << diff << std::endl;
        derivativesGlobal.push_back(diff);
      }
      // end of derivative diagnostics!
      // end of derivative diagnostics!
      // end of derivative diagnostics!
      // end of derivative diagnostics!
      // end of derivative diagnostics!
      // end of derivative diagnostics!


      // avoid outlier hits when applying this cut
      if (!straw_hit._flag.hasAnyProperty(StrawHitFlag::outlier)) {
        if (abs(time_resid) > max_time_res_track) {
          max_time_res_track = abs(time_resid);
        }
      }

      // Convention note:
      // The DoF order is : (planes) dx, dy, dz, a, b, g, followed by (panel) dx, dy, dz, dz,
      // a, b, g This is reflected also in the generated DOCA derivatives.
      std::vector<int> global_dof_labels = generateDOFLabels(straw_id);

      if (_diag > 0 && global_dof_labels.size() != _expected_dofs &&
          derivativesGlobal.size() > _expected_dofs) {
        throw cet::exception("RECO") << "Did not see " << _expected_dofs << " DOF labels"
                                     << " or N of global derivatives was greater than "
                                     << _expected_dofs << " ... Something is wrong!";
      }

      if (derivativesLocal.size() != 5 && use_timeresid) {
        throw cet::exception("RECO")
            << "Did not see 5 local derivatives (corrsp. to 5 fit "
            << "parameters) ... This is weird! (N.B. UseTimeDomain is TRUE)";
      }

      // write the hit to the track buffer
      if (use_timeresid) {
        residuals.emplace_back(time_resid);
        meas_cov(nHits, nHits) = drift_res * drift_res;
      } else {
        residuals.emplace_back(dca_resid);
        meas_cov(nHits, nHits) = drift_res_dca * drift_res_dca;
      }

      for (size_t col = 0; col < 5; ++col) { // FIXME!
        resid_local_derivs(nHits, col) = derivativesLocal[col];
      }

      global_derivs_temp.push_back(
          std::vector<float>(derivativesGlobal.begin(), derivativesGlobal.end()));

      // FIXME!
      local_derivs_temp.push_back(
          std::vector<float>(derivativesLocal.begin(), derivativesLocal.end()));
      labels_temp.push_back(global_dof_labels);

      wrote_hits = true;
      ++nHits;
    }

    if (wrote_hits) {
      // number of hits - 5 track parameters
      ndof = sts._straw_chits.size() - 5;

      if (ndof > 0) {
        pvalue = boost::math::cdf(boost::math::chi_squared(ndof), chisq);
        chisq /= ndof;
      } else {
        chisq = -1;
        pvalue = -1;
        ndof = -1;
        bad_track = true;
      }

      planes_trav = planes_traversed.size();
      panels_trav = panels_traversed.size();

      // track acceptance cuts
      if ((min_plane_traverse != 0 && planes_trav < min_plane_traverse) ||

          (min_panel_traverse_per_plane != 0 &&
           (panels_trav / planes_trav) < min_panel_traverse_per_plane) ||
          (pvalue > max_pvalue) || (max_time_res_track > max_timeres && max_timeres > 0) ||
          (nHits < min_track_hits) || bad_track) {

        if (_diag > 0) {
          std::cout << "track failed quality cuts" << std::endl;
        }

        if (_diag > 2) {
          std::cout << "reason: ";

          if (min_plane_traverse != 0 && planes_trav < min_plane_traverse) {
            std::cout << "not enough planes traversed" << std::endl;
          }

          if (pvalue > max_pvalue) {
            std::cout << "pvalue" << std::endl;
          }

          if (max_time_res_track > max_timeres && max_timeres > 0) {
            std::cout << "max time residual reached (" << max_time_res_track << " > " << max_timeres
                      << ")" << std::endl;
          }

          if (nHits < min_track_hits) {
            std::cout << "hits" << std::endl;
          }

          if (bad_track) {
            std::cout << "bad track" << std::endl;
          }
        }

        continue;
      }
      /*
      Residual Covariance calculation starts here

      FIXME! Still some issues to be ironed out here
      */

      if (_diag > 2) {
        std::cout << "Measurement Covariance:" << std::endl;
        meas_cov.Print();

        std::cout << "Track Covariance:" << std::endl;
        track_cov.Print();
      }

      TMatrixD HC(nHits, 5);
      HC.Mult(resid_local_derivs, track_cov);

      // Using MultT so HCH^T = HC * H^T with H = resid_local_derivs
      TMatrixD HCH(nHits, nHits);
      HCH.MultT(HC, resid_local_derivs);

      meas_cov -= HCH; // V - H C H^T - now holding residual cov

      if (_diag > 2) {
        std::cout << "Residual covariance:" << std::endl;
        meas_cov.Print(); // now holding residual cov
      }
      wrote_hits = false;

      // write hits to buffer
      for (size_t i = 0; i < (size_t)nHits; ++i) {
        residual_err[i] = drift_reso[i]; // sqrt(meas_cov(i, i));
        if (_diag > 2) {
          std::cout << "resid: " << residuals[i] << ", drift res (time): " << residual_err[i]
                    << ", resid cov sqrt(diag): " << sqrt(meas_cov(i, i)) << std::endl;
        }
        if (use_plane_filter && std::find(plane_filter_list.begin(), plane_filter_list.end(),
                                          plane_uid[i]) == plane_filter_list.end()) {
          // we're not interested in this measurement!
        } else {
          millepede->mille(local_derivs_temp[i].size(), local_derivs_temp[i].data(), _expected_dofs,
                           global_derivs_temp[i].data(), labels_temp[i].data(), residuals[i],
                           residual_err[i]);
          wrote_hits = true;
        }

        if (isnan(residual_err[i])) {
          std::cout << "WARNING: sqrt of residual covariance matrix diagonal R_" << i << "," << i
                    << " was NaN! See matrix above. "
                    << "N.B. not using residual covariance." << std::endl;
          // std::cout << "track skipped" << std::endl;
          // bad_track = true;
          // break;
        }
      }

      if (bad_track || !wrote_hits) {
        millepede->kill();
        continue;
      }

      // Write the track buffer to file
      millepede->end();

      if (_diag > 0) {
        diagtree->Fill();
      }

      tracks_written++;
      wrote_track = true;

      if (_diag > 0) {
        std::cout << "wrote track " << tracks_written << std::endl;
      }
    }
  }
  return wrote_track;
}

int AlignTrackCollector::getLabel(int const& object_cls, int const& obj_uid, int const& dof_id) {
  // object class: 0 - 9 - i.e. 1 for planes, 2 for panels
  // object unique id: 0 - 999 supports up to 999 unique objects which is fine for this level of
  // alignment
  // object dof id: 0 - 9

  // 1 000 0
  return object_cls * 10000 + obj_uid * 10 + dof_id;
}

void AlignTrackCollector::analyze(art::Event const& event) {
  StrawResponse const& _srep = srep_h.get(event.id());
  Tracker const& tracker = _proditionsTracker_h.get(event.id());

  auto stH = event.getValidHandle<CosmicTrackSeedCollection>(_costag);

  if (stH.product() == 0) {
    return;
  }

  CosmicTrackSeedCollection const& coscol = *stH.product();
  filter_CosmicTrackSeedCollection(event, tracker, *_tracker, _srep, coscol);
}

}; // namespace mu2e

using mu2e::AlignTrackCollector;
DEFINE_ART_MODULE(AlignTrackCollector);
