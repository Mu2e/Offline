// Ryunosuke O'Neil, 2019
// Module calling upon Mille to set up bootstrap alignment

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
#include "Mu2eUtilities/inc/TwoLinePCA.hh"
#include "RtypesCore.h"
#include "TTree.h"
#include "boost/math/distributions/chi_squared.hpp"
#include "boost/math/distributions/normal.hpp"

#include "art/Framework/Core/EDAnalyzer.h"                 // for EDAnalyzer
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

#include "TAxis.h" // for TAxis
#include "TH1F.h"  // for TH1F

#include "canvas/Utilities/Exception.h" // for Exception
#include "canvas/Utilities/InputTag.h"  // for InputTag

#include "CLHEP/Vector/ThreeVector.h" // for Hep3Vector

#include "cetlib_except/exception.h" // for exception

#include "fhiclcpp/types/Atom.h"                       // for Atom
#include "fhiclcpp/types/Comment.h"                    // for Comment
#include "fhiclcpp/types/Name.h"                       // for Name
#include "fhiclcpp/types/Table.h"                      // for Table::me...
#include "fhiclcpp/types/detail/validationException.h" // for validatio...

#include "Alignment/inc/AlignTrackTypes.hh"    // for AlignTrac...
#include "Alignment/inc/Mille.h"               // for Mille
#include "Alignment/inc/RigidBodyDOCADeriv.hh" // for CosmicTra...

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

    // Histograms for diagnostics
    TH1F* plane_tracks;
    TH1F* plane_residsum;

    TH1F* track_chisq;
    TH1F* track_pvalue;

  public:
    const size_t _dof_per_plane = 6; // dx, dy, dz, a, b, g (translation, rotation)
    const size_t _dof_per_panel = 6; // dx, dy, dz, a, b, g (translation, rotation)
    const size_t _ndof = StrawId::_nplanes * _dof_per_plane + StrawId::_nupanels * _dof_per_panel;

    struct Config {
        using Name = fhicl::Name;
        using Comment = fhicl::Comment;

        fhicl::Atom<int> diaglvl{Name("diagLevel"), Comment("diagnostic level")};

        fhicl::Atom<art::InputTag> costag{Name("CosmicTrackSeedCollection"),
                                          Comment("tag for cosmic track seed collection")};

        fhicl::Atom<std::string> millefile{
            Name("TrackDataOutputFile"), Comment("Output filename for Millepede track data file")};

        fhicl::Atom<std::string> labelsfile{
            Name("LabelsOutputFile"),
            Comment("Output filename for Millepede label ID's (debug only)"), "none"};

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

        fhicl::Atom<double> mindoca{Name("MinDOCA"),
                                    Comment("Require that the drift distance > MinDOCA"), 0.0};


        fhicl::Atom<double> maxtimeres{Name("MaxTimeRes"),
                                    Comment("Require that the maximum time residual on track hits < MaxTimeRes. Setting a negative value does not apply the cut."), -1.0};

        fhicl::Atom<int> mintrackhits{Name("MinTrackSH"),
                                    Comment("Require that the minimum Straw hits in an event > MinTrackSH."),0};


    };
    typedef art::EDAnalyzer::Table<Config> Parameters;

    void beginJob();
    void endJob();
    void beginRun(art::Run const&);
    void analyze(art::Event const&);
    bool filter_CosmicTrackSeedCollection(art::Event const& event, Tracker const& tracker,
                                          StrawResponse const& _srep,
                                          CosmicTrackSeedCollection const& _coscol);
    void writeLabelsFile(Tracker const& aligned_tracker);
    int getLabel(int const&, int const&, int const&);

    AlignTrackCollector(const Parameters& conf)
        : art::EDAnalyzer(conf), _diag(conf().diaglvl()), _costag(conf().costag()),
          _output_filename(conf().millefile()), _labels_filename(conf().labelsfile()),
          track_type(conf().tracktype()), min_plane_traverse(conf().minplanetraverse()),
          min_panel_traverse_per_plane(conf().minpaneltraverse()), max_pvalue(conf().maxpvalue()),
          min_doca(conf().mindoca()), max_timeres(conf().maxtimeres()), min_track_hits(conf().mintrackhits())
    {
        // generate hashtable of plane, or panel number to DOF labels
        // millepede wants arrays of labels

        int nobj = 0;
        // o_cls == 1: planes
        // o_cls == 2: panels
        for (int o_cls = 1; o_cls < 3; o_cls++) // either 1 or 2
        {
            if (o_cls == 1) {
                nobj = StrawId::_nplanes;
            }
            else if (o_cls == 2) {
                nobj = StrawId::_nupanels;
            }
            else
                continue;

            for (int i = 0; i < nobj; i++) {
                std::vector<int> labels;
                for (size_t dof_n = 0; dof_n < 6; dof_n++)
                    labels.push_back(getLabel(o_cls, i, dof_n));

                if (o_cls == 1)
                    plane_dof_labels[i] = std::move(labels);
                else
                    panel_dof_labels[i] = std::move(labels);
            }
        }
        if (_diag > 0) {
            std::cout << "AlignTrackCollector: Total number of plane degrees of freedom = "
                      << StrawId::_nplanes * _dof_per_plane << std::endl;

            std::cout << "AlignTrackCollector: Total number of panel degrees of freedom = "
                      << StrawId::_nupanels * _dof_per_panel << std::endl;
        }

        if (track_type == "CosmicTrackSeedCollection") {
            collect_track = CosmicRecoTrack;
        }
        else {
            throw cet::exception("RECO")
                << "AlignTrackCollector: Cannot collect track type " << track_type << std::endl;
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
    double min_doca;
    double max_timeres;

    int min_track_hits;

    AlignTrackType collect_track;

    std::unique_ptr<Mille> millepede;
    const CosmicTrackSeedCollection* _coscol;
    const Tracker* _tracker;

    size_t tracks_written = 0;

    std::unordered_map<uint16_t, std::vector<int>> plane_dof_labels;
    std::unordered_map<uint16_t, std::vector<int>> panel_dof_labels;

    ProditionsHandle<Tracker> _proditionsTracker_h;
    ProditionsHandle<StrawResponse> srep_h;
};

void AlignTrackCollector::beginJob()
{
    millepede = std::make_unique<Mille>(_output_filename.c_str());

    if (_diag > 0) {
        art::ServiceHandle<art::TFileService> tfs;

        diagtree = tfs->make<TTree>("tracks", "Tracks collected for an alignment iteration");
        diagtree->Branch("nHits", &nHits, "nHits/I");
        diagtree->Branch("doca_resid", &doca_residual, "doca_resid[nHits]/F");
        diagtree->Branch("time_resid", &time_residual, "time_resid[nHits]/F");
        diagtree->Branch("doca_resid_err", &doca_resid_err, "doca_resid_err[nHits]/F");
        diagtree->Branch("drift_res", &drift_reso, "drift_res[nHits]/F");


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

        plane_tracks = tfs->make<TH1F>("plane_trackcount", "Tracks per plane", 36, 0, 36);
        plane_residsum = tfs->make<TH1F>("plane_residusum", "Residual sum per plane", 36, 0, 36);

        track_chisq = tfs->make<TH1F>("track_chisq", "Track ChiSq/NDOF", 30, 0, 4);
        track_pvalue = tfs->make<TH1F>("track_pvalue", "Track p-value", 30, 0, 1);
    }
}

void AlignTrackCollector::writeLabelsFile(Tracker const& aligned_tracker)
{
    // nominal tracking geometry
    Tracker const& nominal_tracker = *GeomHandle<Tracker>();

    // essentially a form of csv output letting you know what labels were
    // assigned to each element
    if (_labels_filename == "none") {
        return;
    }
    std::ofstream label_info_file(_labels_filename);

    label_info_file << "%element_type,uid,dof_idx,global_label" << std::endl;

    for (uint16_t p = 0; p < StrawId::_nplanes; ++p) {
        for (size_t l_idx = 0; l_idx < plane_dof_labels[p].size(); ++l_idx) {
            label_info_file << "Plane," << p << "," << l_idx << "," << plane_dof_labels[p][l_idx]
                            << std::endl;
        }
        label_info_file << "# position diff: "
                        << aligned_tracker.getPlane(p).origin() -
                               nominal_tracker.getPlane(p).origin()
                        << std::endl;


        for (Panel const* panel : aligned_tracker.getPlane(p).getPanels()) {
            uint16_t pa = panel->id().uniquePanel();
            for (size_t l_idx_pa = 0; l_idx_pa < panel_dof_labels[pa].size(); ++l_idx_pa) {
                label_info_file << "Panel," << pa << "," << l_idx_pa << ","
                                << panel_dof_labels[pa][l_idx_pa] << std::endl;
            }
            label_info_file << "# position diff: "
                            << panel->origin() -
                                   nominal_tracker.getPlane(p).getPanel(panel->id()).origin()
                            << std::endl;
            // auto straw0 = panel->getStraw(0);
            // auto straw0pos = straw0.getMidPoint();
            // auto straw0dir = straw0.getDirection();

            // label_info_file << "# generated code Panel Straw0 actual vs. expected position: "
            //     << "actual pos: " << straw0pos << " dir: " << straw0dir
            //     << "expected " << CosmicTrack_DCAalignpos_x(0.5, 0.5, 0,
            //     0,0,0, 0,0,0, 0,0,0,
            //     double plane_x, double plane_y, double plane_z,

            //     double panel_straw0x, double panel_straw0y, double panel_straw0z,

            //     double wire_x, double wire_y, double wire_z,
            //     double wdir_x, double wdir_y, double wdir_z)

        }

        label_info_file << std::endl;
    }

    label_info_file.close();

    std::cout << "AlignTrackCollector: wrote labels to " << _labels_filename << std::endl;
}

void AlignTrackCollector::beginRun(art::Run const& run)
{
    writeLabelsFile(_proditionsTracker_h.get(run.id()));

    return;
}

void AlignTrackCollector::endJob()
{
    // ensure the file is closed once the job finishes
    millepede->~Mille();

    if (_diag > 0) {
        std::cout << "AlignTrackCollector: wrote " << tracks_written << " tracks to "
                  << _output_filename << std::endl;
    }
}

bool AlignTrackCollector::filter_CosmicTrackSeedCollection(art::Event const& event,
                                                           Tracker const& tracker,
                                                           StrawResponse const& _srep,
                                                           CosmicTrackSeedCollection const& coscol)
{
    bool wrote_track = false; // did we write any tracks at all?

    // dedicated to CosmicTrackSeedCollection
    for (CosmicTrackSeed const& sts : coscol) {
        CosmicTrack const& st = sts._track;
        TrkFitFlag const& status = sts._status;

        if (!status.hasAllProperties(TrkFitFlag::helixOK)) {
            continue;
        }

        if (!st.converged || (!st.minuit_converged)) {
            continue;
        }

        if (isnan(st.MinuitParams.A0)) {
            continue;
        }


        std::set<uint16_t> planes_traversed;
        std::set<uint16_t> panels_traversed;

        XYZVec track_pos(st.MinuitParams.A0, 0, st.MinuitParams.B0);
        XYZVec track_dir(st.MinuitParams.A1, -1, st.MinuitParams.B1);

        A0 = st.MinuitParams.A0; // track_pos.X();
        A1 = st.MinuitParams.A1; // track_dir.X();
        B0 = st.MinuitParams.B0;
        B1 = st.MinuitParams.B1;
        T0 = st.MinuitParams.T0;

        GaussianDriftFit fit_object(sts._straw_chits, _srep, &tracker);

        chisq = 0;
        chisq_doca = 0;
        ndof = 0;
        pvalue = 0;
        nHits = 0;

        bool wrote_hits = false; // did we write any hits for the track?
        bool bad_track = false;

        double max_time_res_track = -1;

        // get residuals and their derivatives with respect
        // to all local and global parameters
        // get also plane id hit by straw hits
        for (ComboHit const& straw_hit : sts._straw_chits) {
            // straw and plane info
            StrawId const& straw_id = straw_hit.strawId();
            Straw const& straw = tracker.getStraw(straw_id);
            auto plane_id = straw_id.getPlane();
            auto panel_uuid = straw_id.uniquePanel();
            auto panel_id = straw_id.getPanelId();

            // geometry info
            auto const& plane_origin = tracker.getPlane(plane_id).origin();
            auto const& panel_origin = tracker.getPanel(panel_id).straw0MidPoint();
            auto const& straw_mp = straw.getMidPoint();
            auto const& wire_dir = straw.getDirection().unit();

            // now calculate the derivatives.
            auto derivativesLocal =
                CosmicTrack_DCA_LocalDeriv(A0, B0, A1, B1,
                                           straw_mp.x(), straw_mp.y(), straw_mp.z(),
                                           wire_dir.x(), wire_dir.y(), wire_dir.z(),

                                           plane_origin.x(), plane_origin.y(), plane_origin.z(),
                                           panel_origin.x(), panel_origin.y(), panel_origin.z());

            auto derivativesGlobal =
                CosmicTrack_DCA_GlobalDeriv(A0, B0, A1, B1,
                                            straw_mp.x(), straw_mp.y(), straw_mp.z(),
                                            wire_dir.x(), wire_dir.y(), wire_dir.z(),

                                            plane_origin.x(), plane_origin.y(), plane_origin.z(),
                                            panel_origin.x(), panel_origin.y(), panel_origin.z());

            double resid_tmp = fit_object.DOCAresidual(straw_hit, sts);
            double time_resid = fit_object.TimeResidual(straw_hit, sts);
            double resid_err_tmp = fit_object.DOCAresidualError(straw_hit, sts);


            // FIXME: crude! doesn't belong here!
            CLHEP::Hep3Vector intercept(A0, 0, B0);
            CLHEP::Hep3Vector dir(A1, -1, B1);
            dir = dir.unit();
            TwoLinePCA pca(straw.getMidPoint(), straw.getDirection(), intercept, dir);

            // this is a time
            double drift_res = _srep.driftDistanceError(straw_hit.strawId(), 0, 0, pca.dca());

            chisq += pow(time_resid / drift_res, 2);

            chisq_doca += pow(resid_tmp / resid_err_tmp, 2);

            if (isnan(resid_tmp) || isnan(time_resid) || isnan(drift_res)) {
                bad_track = true;
                continue;
            }

            // FIXME perhaps should get rid of this

            if (pca.dca() < min_doca) {
                continue; // remove hit!
            }

            planes_traversed.insert(plane_id);
            panels_traversed.insert(panel_uuid);

            doca_residual[nHits] = resid_tmp;
            time_residual[nHits] = time_resid;
            doca_resid_err[nHits] = resid_err_tmp;
            pull_doca[nHits] = resid_tmp / resid_err_tmp;
            pull_hittime[nHits] = time_resid / drift_res;

            drift_reso[nHits] = drift_res;

            doca[nHits] = pca.dca();
            time[nHits] = straw_hit.time();
            panel_uid[nHits] = panel_uuid;
            plane_uid[nHits] = plane_id;

            if (_diag > 1) {
                std::cout << "pl" << plane_id << " pa" << panel_uuid << ": resid " << resid_tmp
                          << " +- " << resid_err_tmp << std::endl;
            }

            if (abs(time_resid) > max_time_res_track)
                max_time_res_track = abs(time_resid);

            // FIXME! seems messy!
            std::vector<int> global_dof_labels;
            global_dof_labels.reserve(_dof_per_plane + _dof_per_panel);

            // Convention note:
            // The DoF order is : (planes) dx, dy, dz, a, b, g, followed by (panel) dx, dy, dz, dz,
            // a, b, g This is reflected also in the generated DOCA derivatives.
            global_dof_labels.insert(global_dof_labels.end(), plane_dof_labels[plane_id].begin(),
                                     plane_dof_labels[plane_id].end());
            global_dof_labels.insert(global_dof_labels.end(), panel_dof_labels[panel_uuid].begin(),
                                     panel_dof_labels[panel_uuid].end());

            if (_diag > 0 && global_dof_labels.size() != 12)
                std::cout << "WARNING: we should have 12 labels!" << std::endl;

            // write the hit to the track buffer
            millepede->mille(derivativesLocal.size(), derivativesLocal.data(),
                             derivativesGlobal.size(), derivativesGlobal.data(),
                             global_dof_labels.data(), (float)resid_tmp, (float)resid_err_tmp);

            // diagnostic information
            if (_diag > 0) {
                plane_tracks->Fill(plane_id);
                plane_residsum->Fill(plane_id, resid_tmp);
            }

            wrote_hits = true;
            ++nHits;
        }
        if (wrote_hits) {
            ndof = sts._straw_chits.size() - 5; // 5 track parameters

            if (ndof > 0 && _diag > 0) {
                pvalue = boost::math::cdf(boost::math::chi_squared(ndof), chisq);
                chisq /= ndof;
            }
            else {
                chisq = -1;
                pvalue = -1;
                ndof = -1;
            }

            planes_trav = planes_traversed.size();
            panels_trav = panels_traversed.size();

            // track acceptance cuts
            if ((min_plane_traverse != 0 && planes_trav < min_plane_traverse) ||
                (min_panel_traverse_per_plane != 0 &&
                 (panels_trav / planes_trav) < min_panel_traverse_per_plane) ||
                (pvalue > max_pvalue) ||
                (max_time_res_track > max_timeres && max_timeres > 0) ||
                (nHits < min_track_hits) ||
                bad_track) {
                millepede->kill(); // delete track from buffer

                if (_diag > 0) {
                    std::cout << "track failed quality cuts" << std::endl;
                }

                continue;
            }

            // Write the track buffer to file if we wrote a track
            millepede->end();

            if (_diag > 0)
                diagtree->Fill();
            tracks_written++;
            std::cout << "wrote track " << tracks_written << std::endl;
            wrote_track = true;
        }
    }
    return wrote_track;
}

bool filter_CosmicKalSeedCollection(art::Event const& event, Tracker const& tracker,
                                    StrawResponse const& _srep,
                                    CosmicTrackSeedCollection const& coscol)
{
    return false;
}

int AlignTrackCollector::getLabel(int const& object_cls, int const& obj_uid, int const& dof_id)
{
    // object class: 0 - 9 - i.e. 1 for planes, 2 for panels
    // object unique id: 0 - 999 supports up to 999 unique objects which is fine for this level of
    // alignment
    // object dof id: 0 - 9

    // 1 000 0
    return object_cls * 10000 + obj_uid * 10 + dof_id;
}

void AlignTrackCollector::analyze(art::Event const& event)
{
    StrawResponse const& _srep = srep_h.get(event.id());
    Tracker const& tracker = _proditionsTracker_h.get(event.id());

    switch (collect_track) {

    case CosmicRecoTrack: {
        auto stH = event.getValidHandle<CosmicTrackSeedCollection>(_costag);
        if (stH.product() == 0)
            return;

        CosmicTrackSeedCollection const& coscol = *stH.product();
        filter_CosmicTrackSeedCollection(event, tracker, _srep, coscol);

        break;
    }
    case CosmicKalmanTrack:
        break;
    }
}

}; // namespace mu2e

using mu2e::AlignTrackCollector;
DEFINE_ART_MODULE(AlignTrackCollector);
