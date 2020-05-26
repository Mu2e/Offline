// Ryunosuke O'Neil, 2019
// roneil@fnal.gov
// ryunoneil@gmail.com

// A module to SELECT Cosmic NoField events that yield a track and pass some cuts.

// Consult README.md for more information

#include <algorithm> // for max, all_of
#include <cmath>     // for isnan
#include <cstddef>   // for size_t
#include <exception>
#include <iostream>      // for operator<<
#include <memory>        // for unique_ptr
#include <stdint.h>      // for uint16_t
#include <string>        // for string
#include <unordered_map> // for unordered...
#include <utility>       // for move
#include <vector>        // for vector<>:...

#include "CosmicReco/inc/PDFFit.hh"
#include "GeneralUtilities/inc/BitMap.hh"
#include "GeometryService/inc/GeomHandle.hh"
#include "Minuit2/MnUserCovariance.h"
#include "Mu2eUtilities/inc/TwoLinePCA.hh"

#include "boost/math/distributions/chi_squared.hpp"
#include "boost/math/distributions/normal.hpp"

#include "art/Framework/Core/EDFilter.h"                   // for EDAnalyzer
#include "art/Framework/Core/ModuleMacros.h"               // for DEFINE_AR...
#include "art/Framework/Core/ProducerTable.h"              // for ProducerT...
#include "art/Framework/Principal/Event.h"                 // for Event
#include "art/Framework/Principal/Handle.h"                // for ValidHandle
#include "art/Framework/Services/Registry/ServiceHandle.h" // for ServiceHa...
#include "art_root_io/TFileService.h"                      // for TFileService

#include "TrackerConditions/inc/StrawResponse.hh" // for StrawResp...
#include "TrackerGeom/inc/Panel.hh"               // for Panel
#include "TrackerGeom/inc/Plane.hh"               // for Plane
#include "TrackerGeom/inc/Straw.hh"               // for Straw
#include "TrackerGeom/inc/Tracker.hh"             // for Tracker

#include "ProditionsService/inc/ProditionsHandle.hh" // for Prodition...
#include "RecoDataProducts/inc/ComboHit.hh"          // for ComboHit
#include "RecoDataProducts/inc/CosmicTrack.hh"       // for CosmicTrack
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"   // for CosmicTra...
#include "RecoDataProducts/inc/TrkFitFlag.hh"        // for TrkFitFlag

#include "DataProducts/inc/StrawId.hh" // for StrawId
#include "DataProducts/inc/XYZVec.hh"  // for toXYZVec

#include "Mu2eUtilities/inc/TwoLinePCA_XYZ.hh" // for TwoLinePC...

#include "RtypesCore.h"
#include "TAxis.h" // for TAxis
#include "TH1F.h"  // for TH1F
#include "TMatrixDSym.h"
#include "TTree.h"

#include "canvas/Utilities/Exception.h" // for Exception
#include "canvas/Utilities/InputTag.h"  // for InputTag

#include "CLHEP/Vector/ThreeVector.h" // for Hep3Vector

#include "cetlib_except/exception.h" // for exception

#include "fhiclcpp/types/Atom.h"                       // for Atom
#include "fhiclcpp/types/Comment.h"                    // for Comment
#include "fhiclcpp/types/Name.h"                       // for Name
#include "fhiclcpp/types/Table.h"                      // for Table::me...
#include "fhiclcpp/types/detail/validationException.h" // for validatio...

#include "TrackerAlignment/inc/AlignmentDerivatives.hh" // for CosmicTra...
#include "TrackerAlignment/inc/Mille.h"                 // for Mille

namespace art {
class Run;
} // namespace art

using namespace mu2e;
namespace mu2e {

class AlignTrackSelector : public art::EDFilter {
private:
  static constexpr int MAX_NHITS = 100;

  Int_t nHits;

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

    fhicl::Atom<double> mindoca{Name("MinDOCA"),
                                Comment("Require that the drift distance > MinDOCA"), 0.0};

    fhicl::Atom<double> maxtimeres{Name("MaxTimeRes"),
                                   Comment("Require that the maximum ABSOLUTE time residual over "
                                           "all track hits < MaxTimeRes. Setting a "
                                           "negative value does not apply the cut."),
                                   -1.0};

    fhicl::Atom<int> mintrackhits{
        Name("MinTrackSH"), Comment("Require that the minimum straw hits in a track > MinTrackSH."),
        0};
  };
  typedef art::EDFilter::Table<Config> Parameters;

  void beginJob();
  void endJob();
  bool beginRun(art::Run&);
  bool filter(art::Event&);
  bool filter_CosmicTrackSeedCollection(art::Event const& event, Tracker const& tracker,
                                        StrawResponse const& _srep,
                                        CosmicTrackSeedCollection const& _coscol);

  int getLabel(int const&, int const&, int const&);
  std::vector<int> generateDOFLabels(StrawId const& strw);

  AlignTrackSelector(Parameters const& conf) :
      art::EDFilter(conf), _diag(conf().diaglvl()), _costag(conf().costag()),

      min_plane_traverse(conf().minplanetraverse()),
      min_panel_traverse_per_plane(conf().minpaneltraverse()), max_pvalue(conf().maxpvalue()),
      min_doca(conf().mindoca()), max_timeres(conf().maxtimeres()),
      min_track_hits(conf().mintrackhits()) {

    consumes<CosmicTrackSeedCollection>(_costag);
  }

  virtual ~AlignTrackSelector() {}

  Config _conf;

  int _diag;
  art::InputTag _costag;

  int min_plane_traverse;
  int min_panel_traverse_per_plane;
  double max_pvalue;
  double min_doca;
  double max_timeres;
  int min_track_hits;

  const CosmicTrackSeedCollection* _coscol;
  const Tracker* _tracker;

  size_t tracks_written = 0;

  ProditionsHandle<Tracker> _proditionsTracker_h;
  ProditionsHandle<StrawResponse> srep_h;
};

void AlignTrackSelector::beginJob() {}

bool AlignTrackSelector::beginRun(art::Run& run) { return true; }

void AlignTrackSelector::endJob() {
  if (_diag > 0) {
    std::cout << "AlignTrackSelector: selected " << tracks_written << " tracks " << std::endl;
  }
}

bool AlignTrackSelector::filter_CosmicTrackSeedCollection(art::Event const& event,
                                                          Tracker const& tracker,
                                                          StrawResponse const& _srep,
                                                          CosmicTrackSeedCollection const& coscol) {
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

    XYZVec track_pos(st.MinuitParams.A0, 0, st.MinuitParams.B0);
    XYZVec track_dir(st.MinuitParams.A1, -1, st.MinuitParams.B1);

    A0 = st.MinuitParams.A0;
    A1 = st.MinuitParams.A1;
    B0 = st.MinuitParams.B0;
    B1 = st.MinuitParams.B1;
    T0 = st.MinuitParams.T0;

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

      auto plane_id = straw_id.getPlane();
      auto panel_uuid = straw_id.uniquePanel();

      double resid_tmp = fit_object.DOCAresidual(straw_hit, sts);
      double time_resid = fit_object.TimeResidual(straw_hit, sts);
      
      CLHEP::Hep3Vector intercept(A0, 0, B0);
      CLHEP::Hep3Vector dir(A1, -1, B1);
      dir = dir.unit();
      TwoLinePCA pca(straw.getMidPoint(), straw.getDirection(), intercept, dir);

      // FIXME! this is a time, not distance
      double drift_res = _srep.driftDistanceError(straw_hit.strawId(), 0, 0, pca.dca());
      double resid_err_tmp = _srep.driftTimeToDistance(
          straw_hit.strawId(), drift_res, 0); // fit_object.DOCAresidualError(straw_hit, sts);

      // FIXME! use newly implemented chisq function in fit object
      chisq += pow(time_resid / drift_res, 2);
      chisq_doca += pow(resid_tmp / resid_err_tmp, 2);

      if (isnan(resid_tmp) || isnan(time_resid) || isnan(drift_res)) {
        bad_track = true;
        continue;
      }

      planes_traversed.insert(plane_id);
      panels_traversed.insert(panel_uuid);

      if (_diag > 2) {
        std::cout << "pl" << plane_id << " pa" << panel_uuid << ": resid " << resid_tmp << " +- "
                  << resid_err_tmp << std::endl;
      }

      // avoid outlier hits when applying this cut
      if (!straw_hit._flag.hasAnyProperty(StrawHitFlag::outlier)) {
        if (abs(time_resid) > max_time_res_track) {
          max_time_res_track = abs(time_resid);
        }
      }

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

          (nHits < min_track_hits) ||

          bad_track) {

        if (_diag > 0) {
          std::cout << "track failed quality cuts" << std::endl;
        }
        bad_track = true;
        continue;
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

bool AlignTrackSelector::filter(art::Event& event) {
  StrawResponse const& _srep = srep_h.get(event.id());
  Tracker const& tracker = _proditionsTracker_h.get(event.id());

  auto stH = event.getValidHandle<CosmicTrackSeedCollection>(_costag);

  if (stH.product() == 0) {
    return false;
  }
  CosmicTrackSeedCollection const& coscol = *stH.product();

  return filter_CosmicTrackSeedCollection(event, tracker, _srep, coscol);
}

}; // namespace mu2e

using mu2e::AlignTrackSelector;
DEFINE_ART_MODULE(AlignTrackSelector);
