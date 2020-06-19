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

#include "CosmicReco/inc/MinuitDriftFitter.hh"
#include "CosmicReco/inc/PDFFit.hh"
#include "DbTables/inc/TrkAlignPanel.hh"
#include "DbTables/inc/TrkAlignPlane.hh"
#include "GeneralUtilities/inc/BitMap.hh"
#include "GeneralUtilities/inc/HepTransform.hh"
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
#include "TrackerAlignment/inc/MilleDataWriter.hh"
#include "TrackerAlignment/inc/AlignmentUtilities.hh"

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

    fhicl::Atom<std::string> constrainstrategy{
        Name("ConstrainStrategy"),
        Comment("Which constraint strategy to use. Either 'Fix' or 'Average'"), "Fix"};
    
    fhicl::Sequence<int> fixplane{
        Name("FixPlane"),
        Comment("Planes to fix. The parameters are set to the proditions value."),
    };

    fhicl::Atom<std::string> millefile{Name("MilleFile"),
                                       Comment("Output filename for Millepede track data file")};

    fhicl::Atom<bool> millefilegzip{Name("GzipCompression"),
                                    Comment("Enable gzip compression for millepede output file"),
                                    false};

    fhicl::Atom<std::string> steerfile{
        Name("SteerFile"), Comment("Output filename for Millepede steering file"), "steer.txt"};

    fhicl::Atom<std::string> paramfile{
        Name("ParamFile"), Comment("Output filename for Millepede parameters file"), "params.txt"};

    fhicl::Atom<std::string> constrfile{Name("ConstrFile"),
                                        Comment("Output filename for Millepede constraints file"),
                                        "constr.txt"};

    fhicl::Atom<bool> usenumerical{
        Name("UseNumericalDiffn"),
        Comment("Whether or not to use numerical derivatives. Default is false."), false};

    fhicl::Sequence<std::string> mpsteers{
        Name("SteeringOpts"),
        Comment("Additional configuration options to add to generated steering file."),
    };

    fhicl::Atom<double> errorscale{
        Name("ErrorScale"),
        Comment("Scale all hit measurement errors by ErrorScale."), 1.0
    };

    fhicl::Atom<bool> enableCV{
        Name("EnableLOOCVFitting"),
        Comment("Whether to enable LOOCV track fitting to obtain 'unbiased' residuals. Default is true."), true};

  };

  typedef art::EDAnalyzer::Table<Config> Parameters;

  Config _conf;

  int _diag;
  art::InputTag _costag;
  std::string _labels_filename;

  int min_plane_traverse;
  int min_panel_traverse_per_plane;
  double max_pvalue;
  double max_timeres;
  int min_track_hits;
  bool use_timeresid;
  bool no_panel_dofs;
  bool no_plane_rotations;
  bool use_db;

  bool gzip_compress;
  std::string mille_filename;
  std::string steer_filename;
  std::string param_filename;
  std::string constr_filename;

  bool use_numeric_derivs;
  bool wroteMillepedeParams;
  MilleDataWriter<double> mille_file;
  std::vector<std::string> steer_lines;

  double error_scale;
  bool use_unbiased_res;


  std::string constrain_strat;
  std::vector<int> fixed_planes;
  

  const CosmicTrackSeedCollection* _coscol;
  const Tracker* _tracker;

  size_t tracks_written = 0;

  ProditionsHandle<Tracker> _proditionsTracker_h;
  ProditionsHandle<StrawResponse> srep_h;

  std::unique_ptr<DbHandle<TrkAlignPlane>> _trkAlignPlane_h;
  std::unique_ptr<DbHandle<TrkAlignPanel>> _trkAlignPanel_h;

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

  void writeMillepedeSteering();
  void writeMillepedeConstraints();
  void writeMillepedeParams(TrkAlignPlane const& alignConstPlanes,
                            TrkAlignPanel const& alignConstPanels);
  bool isDOFenabled(int object_class, int object_id, int dof_n);
  bool isDOFfixed(int object_class, int object_id, int dof_n);
  
double CosmicTrack_RealDCA(
  StrawId const& strawId,

  double const& a0, double const& b0, 
  double const& a1, double const& b1, 
  double const& t0, 

  double const& plane_dx, double const& plane_dy, double const& plane_dz, 
  double const& plane_a, double const& plane_b, double const& plane_g, 

  double const& panel_dx, double const& panel_dy, double const& panel_dz, 
  double const& panel_a, double const& panel_b, double const& panel_g, 

  StrawResponse const& _srep, int ambiguity);


  virtual ~AlignTrackCollector() {}

  AlignTrackCollector(const Parameters& conf) :
      art::EDAnalyzer(conf), 
      _diag(conf().diaglvl()),
      _costag(conf().costag()),
      min_plane_traverse(conf().minplanetraverse()),
      min_panel_traverse_per_plane(conf().minpaneltraverse()), 
      max_pvalue(conf().maxpvalue()),
      max_timeres(conf().maxtimeres()), 
      min_track_hits(conf().mintrackhits()),
      use_timeresid(conf().usetimeresid()), 
      no_panel_dofs(conf().nopaneldofs()),
      no_plane_rotations(conf().noplanerotations()), 

      gzip_compress(conf().millefilegzip()), 
      mille_filename(conf().millefile()),
      steer_filename(conf().steerfile()), 
      param_filename(conf().paramfile()),
      constr_filename(conf().constrfile()),

      use_numeric_derivs(conf().usenumerical()),

      wroteMillepedeParams(false), 
      mille_file(mille_filename, gzip_compress),
      steer_lines(conf().mpsteers()),
      error_scale(conf().errorscale()),
      use_unbiased_res(conf().enableCV()),
      constrain_strat(conf().constrainstrategy()),
      fixed_planes(conf().fixplane()) {

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

      std::cout << "AlignTrackCollector: compression is "
                << (gzip_compress ? "enabled" : "disabled") << std::endl;
    }
  }
};

void AlignTrackCollector::beginJob() {
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

void AlignTrackCollector::writeMillepedeSteering() {
  std::ofstream output_file(steer_filename);

  output_file << "! Steering file generated by AlignTrackCollector" << std::endl
              << "Cfiles" << std::endl
              << param_filename << std::endl
              << constr_filename << std::endl
              << mille_filename << std::endl
              << std::endl;

  for (std::string const& line : steer_lines) {
    output_file << line << std::endl;
  }
  
  output_file << std::endl
              << "method inversion 10 0.001" << std::endl
              << "end" << std::endl;

}

bool AlignTrackCollector::isDOFfixed(int object_class, int object_id, int dof_n) {
  if (object_class != 1) { return false; } // FIXME!
  return std::find(fixed_planes.begin(), fixed_planes.end(), object_id) != fixed_planes.end();
}

bool AlignTrackCollector::isDOFenabled(int object_class, int object_id, int dof_n) {
  if (object_class == 2 && no_panel_dofs) {
    return false;
  }
  if (object_class == 1 && no_plane_rotations && dof_n > 2) {
    return false;
  }
  return true;
}

void AlignTrackCollector::writeMillepedeConstraints() {


  std::ofstream output_file(constr_filename);
  output_file << "! Generated by AlignTrackCollector" << std::endl;
  output_file << "! Constraint strategy: " << constrain_strat << std::endl;

  if (constrain_strat != "Average") { // this function is called once, so the string comp happens once
    return;
  }

  output_file << "! plane constraints:" << std::endl;

  // plane translations, rotations
  for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
    output_file << "Constraint    0" << std::endl;

    for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {

      if (!isDOFenabled(1, p, dof_n)) {
        continue;
      }

      // label(int)   coefficient of DOF corresp. to label in the sum
      output_file << getLabel(1, p, dof_n) << "    1" << std::endl;
    }
  }

  output_file << "! panels follow" << std::endl;

  // panel translations, rotations
  for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
    output_file << "Constraint    0" << std::endl;

    for (uint16_t p = 0; p < StrawId::_nupanels; ++p) {
      if (!isDOFenabled(1, p, dof_n)) {
        continue;
      }

      // label(int)   coefficient of DOF corresp. to label in the sum
      output_file << getLabel(2, p, dof_n) << "    " << "  1" << std::endl;

    }
    output_file << std::endl;
  }


}

void AlignTrackCollector::writeMillepedeParams(TrkAlignPlane const& alignConstPlanes,
                                               TrkAlignPanel const& alignConstPanels) {
  // write a params.txt telling millepede what the alignment constants
  // were for this track collection iteration

  std::ofstream output_file(param_filename);
  output_file << "! Generated by AlignTrackCollector" << std::endl;
  output_file << "! Constraint strategy: " << constrain_strat << std::endl;
  output_file << "! columns: label, alignment parameter start value, "
                 "presigma (-ve: fixed, 0: variable)" << std::endl;

  output_file << "Parameter" << std::endl;
  output_file << "! Plane DOFs:" << std::endl;

  for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
    auto const& rowpl = alignConstPlanes.rowAt(p);

    std::vector<float> pl_consts{rowpl.dx(), rowpl.dy(), rowpl.dz(),
                                 rowpl.rx(), rowpl.ry(), rowpl.rz()};

    for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
      if (!isDOFenabled(1, p, dof_n)) {
        continue;
      }

      // label(int)   initial value (float)   presigma (-ve: fixed, 0: variable)
      output_file << getLabel(1, p, dof_n) << "  " << pl_consts[dof_n] << "  " 
                  << (isDOFfixed(1, p, dof_n) ? "-1" : "0") << std::endl;
    }
  }
  output_file << std::endl << "! end plane DOFs" << std::endl;

  output_file << std::endl << "! Panel DOFs:" << std::endl;

  for (uint16_t p = 0; p < StrawId::_nupanels; ++p) {
    auto const& rowpa = alignConstPanels.rowAt(p);

    std::vector<float> pa_consts{rowpa.dx(), rowpa.dy(), rowpa.dz(),
                                 rowpa.rx(), rowpa.ry(), rowpa.rz()};

    for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
      if (!isDOFenabled(2, p, dof_n)) {
        continue;
      }

      // label(int)   initial value (float)   presigma (-ve: fixed, 0: variable)
      output_file << getLabel(2, p, dof_n) << "  " << pa_consts[dof_n] << "  " << (isDOFfixed(2, p, dof_n) ? "-1" : "0") << std::endl;
    }
  }
  output_file << std::endl << "! end panel DOFs" << std::endl;

  output_file.close();
}

int AlignTrackCollector::getLabel(int const& object_cls, int const& obj_uid, int const& dof_id) {
  // object class: 0 - 9 - i.e. 1 for planes, 2 for panels
  // object unique id: 0 - 999 supports up to 999 unique objects which is fine for this level of
  // alignment
  // object dof id: 0 - 9

  // (future) a straw label might be:
  // id  plane panel straw dof
  // 3   00    0     00    0
  // min
  // 3 00 0 00 0
  // max
  // 3,356,959

  // 1 000 0

  return object_cls * 10000 + obj_uid * 10 + dof_id;
}


std::vector<int> AlignTrackCollector::generateDOFLabels(StrawId const& strw) {
  return generateDOFLabels(strw.getPlane(), strw.uniquePanel());
}

std::vector<int> AlignTrackCollector::generateDOFLabels(uint16_t plane, uint16_t panel) {
  std::vector<int> labels;
  labels.reserve(_dof_per_plane + _dof_per_panel);

  for (size_t dof_n = 0; dof_n < _dof_per_plane; dof_n++) {
    if (!isDOFenabled(1, plane, dof_n)) {
      continue;
    }
    labels.push_back(getLabel(1, plane, dof_n));
  }
  for (size_t dof_n = 0; dof_n < _dof_per_panel; dof_n++) {
    if (!isDOFenabled(2, panel, dof_n)) {
      continue;
    }
    labels.push_back(getLabel(2, panel, dof_n));
  }
  
  return labels;
}

void AlignTrackCollector::beginRun(art::Run const& run) {
  GeomHandle<Tracker> track;
  _tracker = track.get();
}

void AlignTrackCollector::endJob() {

  if (_diag > 0) {
    std::cout << "AlignTrackCollector: wrote " 
              << tracks_written << " tracks to " 
              << mille_filename
              << std::endl;
  }

  writeMillepedeConstraints();
  writeMillepedeSteering();
}

bool AlignTrackCollector::filter_CosmicTrackSeedCollection(
    art::Event const& event, Tracker const& alignedTracker, Tracker const& nominalTracker,
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

    if (!sts._status.hasAllProperties(TrkFitFlag::helixOK)) {
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

    AlignmentUtilities::CosmicTimeTrack original_track { A0, B0, A1, B1, T0 };

    // used/re-used for unbiased residual fits
    AlignmentUtilities::CosmicTimeTrack track { A0, B0, A1, B1, T0 };

    GaussianDriftFit fit_object(sts._straw_chits, _srep, &alignedTracker);

    chisq = 0;
    chisq_doca = 0;
    ndof = 0;
    pvalue = 0;
    nHits = 0;

    // for the max timeresidual track quality cut
    double max_time_res_track = -1;

    bool wrote_hits = false;
    bool bad_track = false;

    std::vector<double> residuals;
    std::vector<std::vector<double>> global_derivs_temp;
    std::vector<std::vector<double>> local_derivs_temp;
    std::vector<std::vector<int>> labels_temp;


    // get residuals and their derivatives with respect
    // to all local and global parameters

    for (ComboHit const& straw_hit : sts._straw_chits) {
      // straw and plane info
      StrawId const& straw_id = straw_hit.strawId();
      Straw const& alignedStraw = alignedTracker.getStraw(straw_id);

      auto plane_id = straw_id.getPlane();
      auto panel_id = straw_id.uniquePanel();

      if (use_unbiased_res)
      {
        // re-fit without this hit.
        fit_object.setExcludeHit(nHits);

        std::vector<double> pars(original_track.as_vector());
        std::vector<double> errors {
          sts._track.MinuitParams.deltaA0,
          sts._track.MinuitParams.deltaB0,
          sts._track.MinuitParams.deltaA1,
          sts._track.MinuitParams.deltaB1,
          sts._track.MinuitParams.deltaT0
        };
        std::vector<double> cov;
        bool converged;

        MinuitDriftFitter::DoDriftTimeFit(
          pars, errors, cov, converged, fit_object);

        // update track variable
        if (converged) {
          track.setParams(pars);
        } else {
          std::cerr << "WARNING: unconverged track for this hit. Using original track ..." << std::endl;
          track.setParams(original_track);
        }
      }

      // hit DOCA, DOCA/TOCA residuals, errors
      TwoLinePCA pca(track.intercept(), track.direction(), 
          alignedStraw.getMidPoint(), alignedStraw.getDirection());

      double driftvel = _srep.driftInstantSpeed(straw_id, pca.dca(), 0);
      double dca_resid = fit_object.DOCAresidual(straw_hit, track.as_vector());
      double drift_res_dca = error_scale * _srep.driftDistanceError(straw_hit.strawId(), 0, 0, pca.dca());

      double time_resid = fit_object.TimeResidual(straw_hit, track.as_vector());
      double drift_res = error_scale * _srep.driftTimeError(straw_hit.strawId(), 0, 0, pca.dca());

      // alignment constants for derivatives
      auto const& rowpl = alignConsts_planes.rowAt(plane_id);
      auto const& rowpa = alignConsts_panels.rowAt(panel_id);

      // now calculate the derivatives
      std::vector<double> derivativesLocal, derivativesGlobal;

      if (use_numeric_derivs) {
        std::tie(derivativesLocal, derivativesGlobal) = AlignmentUtilities::numericalDerivatives(
          track, straw_id, rowpl, rowpa, nominalTracker, _srep);
      }
      else {
        std::tie(derivativesLocal, derivativesGlobal) = AlignmentUtilities::analyticalDerivatives(
          track, straw_id, rowpl, rowpa, nominalTracker, driftvel);
      }

      chisq += pow(time_resid / drift_res, 2);
      chisq_doca += pow(dca_resid / drift_res_dca, 2);

      if (isnan(dca_resid) || isnan(time_resid) || isnan(drift_res)) {
        bad_track = true;
        continue;
      }

      planes_traversed.insert(plane_id);
      panels_traversed.insert(panel_id);

      // fill tree info
      doca_residual[nHits] = dca_resid;
      time_residual[nHits] = time_resid;
      doca_resid_err[nHits] = drift_res_dca;
      pull_doca[nHits] = dca_resid / drift_res_dca;
      pull_hittime[nHits] = time_resid / drift_res;
      drift_reso[nHits] = drift_res;
      doca[nHits] = pca.dca();
      time[nHits] = straw_hit.time();
      panel_uid[nHits] = panel_id;
      plane_uid[nHits] = plane_id;

      if (_diag > 3) {
        std::cout << "pl" << plane_id << " pa" << panel_id 
                  << ": dcaresid " << dca_resid << " +- "
                  << drift_res_dca << std::endl
                  << "timeresid " << time_resid << " +- " 
                  << drift_res << std::endl;
      }

      // Time residual cut
      // avoid outlier hits when applying this cut
      if (!straw_hit._flag.hasAnyProperty(StrawHitFlag::outlier)) {
        if (abs(time_resid) > max_time_res_track) {
          max_time_res_track = abs(time_resid);
        }
      }

      if (_diag > 4) {
        AlignmentUtilities::diagPrintHit(
          track, time_resid, drift_res, derivativesLocal, derivativesGlobal, straw_id);

        if (!AlignmentUtilities::testDerivatives( /// FIXME!
              pca, alignedTracker, track, straw_id, rowpl, rowpa, nominalTracker, _srep)) {
          std::cout << "----------------------------------" << std::endl;
          std::cout
              << "WARNING! AlignmentDerivatives are inconsistent! Please validate."
              << std::endl;
          std::cout << "----------------------------------" << std::endl;

          // throw cet::exception("ALIGNMENT") << "Output of generated functions"
          //   << "(AlignmentDerivatives) are not consistent within expected tolerance! Please validate.";
        }
      }
      
      std::vector<int> labels = generateDOFLabels(straw_id);

      if (_diag > 0 && labels.size() != derivativesGlobal.size() &&
          derivativesGlobal.size() < _expected_dofs) {
        throw cet::exception("RECO") << "N global derivatives != N labels"
                                     << " or N of global derivatives was less than "
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
      } else {
        residuals.emplace_back(dca_resid);
      }

      derivativesGlobal.resize(_dof_per_plane + _dof_per_panel);
      global_derivs_temp.emplace_back(derivativesGlobal);
      local_derivs_temp.emplace_back(derivativesLocal);
      labels_temp.push_back(labels);

      wrote_hits = true;
      ++nHits;
    }

    if (wrote_hits) {
      // number of hits - 5 track parameters
      ndof = sts._straw_chits.size() - 5;

      if (ndof > 0) {
        pvalue = 1.0-boost::math::cdf(boost::math::chi_squared(ndof), chisq);
        chisq /= ndof;
      } else {
        bad_track = true;
      }

      planes_trav = planes_traversed.size();
      panels_trav = panels_traversed.size();

      // track acceptance cuts
      if ((min_plane_traverse != 0 && 
            planes_trav < min_plane_traverse) ||
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
        std::cout << std::endl;
        continue;
      }

      // write hits to buffer
      for (size_t i = 0; i < (size_t)nHits; ++i) {
        residual_err[i] = drift_reso[i];
        mille_file.pushHit(local_derivs_temp[i], global_derivs_temp[i], labels_temp[i],
                            residuals[i], residual_err[i]);
      }

      if (bad_track || !wrote_hits) {
        if (_diag > 0) {
          std::cout << "killed track " << std::endl;
          std::cout << std::endl;
        }
        mille_file.clear();
        continue;
      }

      // Write the track buffer to file
      mille_file.flushTrack();

      if (_diag > 0) {
        diagtree->Fill();
      }

      tracks_written++;
      wrote_track = true;

      if (_diag > 0) {
        std::cout << "wrote track " << tracks_written << std::endl;
        std::cout << std::endl;
      }
    }
  }
  return wrote_track;
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
