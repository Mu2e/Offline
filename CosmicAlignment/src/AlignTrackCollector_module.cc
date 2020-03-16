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

#include "boost/math/distributions/chi_squared.hpp"
#include "boost/math/distributions/normal.hpp"

#include "art/Framework/Core/EDAnalyzer.h"                   // for EDAnalyzer
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

#include "CosmicAlignment/inc/AlignTrackTypes.hh"    // for AlignTrac...
#include "CosmicAlignment/inc/Mille.h"               // for Mille
#include "CosmicAlignment/inc/RigidBodyDOCADeriv.hh" // for CosmicTra...

namespace art
{
    class Run;
} // namespace art

using namespace mu2e;
namespace mu2e
{

    class AlignTrackCollector : public art::EDAnalyzer
    {
        private:
            // Histograms for diagnostic purposes
            TH1F* residuum;

            TH1F* plane_tracks;
            TH1F* plane_residsum;
            TH1F* plane_seedtracks;
            TH1F* plane_seedresidualsum;

            TH1F* track_chisq;
            TH1F* track_pvalue;

            TH1F* seedtrack_chisq;
            TH1F* seedtrack_pvalue;

            TH1F* resid_err;
            TH1F* doca_h;

            TH1F* drift_dist;
            TH1F* drift_time;

        public:
            const size_t _dof_per_plane = 6; // dx, dy, dz, a, b, g (translation, rotation)
            const size_t _dof_per_panel = 6; // dx, dy, dz, a, b, g (translation, rotation)
            const size_t _ndof = StrawId::_nplanes * _dof_per_plane + StrawId::_nupanels * _dof_per_panel;

            struct Config
            {
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

                fhicl::Atom<bool> seedonly {
                    Name("SeedOnly"),
                        Comment("Use seed tracks only"), false
                };
            };
            typedef art::EDAnalyzer::Table<Config> Parameters;

            void beginJob();
            void endJob();
            void beginRun(art::Run const&);
            void analyze(art::Event const&);
            bool filter_CosmicTrackSeedCollection(art::Event const& event, Tracker const& tracker,
                    StrawResponse const& _srep,
                    CosmicTrackSeedCollection const& _coscol);
            void writeLabelsFile();
            int getLabel(int const&, int const&, int const&);

            AlignTrackCollector(const Parameters& conf)
                : art::EDAnalyzer(conf), _diag(conf().diaglvl()), _costag(conf().costag()),
                _output_filename(conf().millefile()), _labels_filename(conf().labelsfile()),
                track_type(conf().tracktype()), seed_only(conf().seedonly())
        {
            // generate hashtable of plane, or panel number to DOF labels
            // we prepare them like this because millepede wants arrays of labels
            // (I'd rather cache them in memory than re-calculate them for every single track)
            int nobj = 0;
            for (int o_cls = 1; o_cls < 3; o_cls++) // labels start at 1.
            {
                if (o_cls == 1)
                {
                    nobj = StrawId::_nplanes;
                }
                else if (o_cls == 2)
                {
                    nobj = StrawId::_nupanels;
                }
                else continue;
                for (uint16_t i = 0; i < nobj; i++)
                {
                    std::vector<int> labels;
                    for (size_t dof_n = 0; dof_n < 6; dof_n++)
                        labels.push_back(getLabel(o_cls, i, dof_n));

                    if (o_cls == 1)
                        plane_dof_labels[i] = std::move(labels);
                    else
                        panel_dof_labels[i] = std::move(labels);
                }
            }
            if (_diag > 0)
            {
                std::cout << "AlignTrackCollector: Total number of plane degrees of freedom = "
                    << StrawId::_nplanes * _dof_per_plane << std::endl;

                std::cout << "AlignTrackCollector: Total number of panel degrees of freedom = "
                    << StrawId::_nupanels * _dof_per_panel << std::endl;
            }

            if (track_type == "CosmicTrackSeedCollection")
            {
                collect_track = CosmicRecoTrack;
            }
            else
            {
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
            bool seed_only;

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

        writeLabelsFile();

        if (_diag > 0)
        {
            // TODO: extend these diagnostics  (could fill a TTree?)
            art::ServiceHandle<art::TFileService> tfs;
            residuum = tfs->make<TH1F>("residuum", "Straw Hit Residuals ", 100, -40, 25);
            residuum->GetXaxis()->SetTitle("Residual (DOCA - Estimated Drift Distance) (mm)");

            resid_err = tfs->make<TH1F>("resid_err", "Straw Hit Drift distance error", 100, -5, 5);

            plane_seedtracks = tfs->make<TH1F>("plane_seedtrackcount", "Tracks per plane", 36, 0, 36);
            // plane_seedresidualsum = tfs->make<TH1F>("plane_seedresidsum", "Residual sum per plane",
            // 36, 0, 36);
            plane_tracks = tfs->make<TH1F>("plane_trackcount", "Tracks per plane", 36, 0, 36);
            plane_residsum = tfs->make<TH1F>("plane_residusum", "Residual sum per plane", 36, 0, 36);

            // seedtrack_chisq = tfs->make<TH1F>("seedtrack_chisq", "Seed Track chi squared/ndof", 30,
            // 0, 4);  seedtrack_pvalue = tfs->make<TH1F>("seedtrack_pvalue", "Seed Track P Value", 30, 0,
            // 1);
            drift_time = tfs->make<TH1F>("drift_time", "Drift Time", 100, 0, 300); // ns
            drift_dist = tfs->make<TH1F>("drift_dist", "Drift Dist", 50, 0, 50); // mm
            doca_h = tfs->make<TH1F>("doca", "DOCA", 50, -4, 4);

            track_chisq = tfs->make<TH1F>("track_chisq", "Track ChiSq/NDOF", 30, 0, 4);
            track_pvalue = tfs->make<TH1F>("track_pvalue", "Track p-value", 30, 0, 1);
        }
    }

    void AlignTrackCollector::writeLabelsFile()
    {
        // essentially a form of csv output letting you know what labels were
        // assigned to each element
        if (_labels_filename == "none")
            return;

        std::ofstream label_info_file(_labels_filename);

        label_info_file << "%element_type,uid,dof_idx,global_label" << std::endl;

        for (uint16_t p = 0; p < StrawId::_nplanes; ++p)
            for (size_t l_idx = 0; l_idx < plane_dof_labels[p].size(); ++l_idx)
                label_info_file << "Plane," << p << "," << l_idx << "," << plane_dof_labels[p][l_idx]
                    << std::endl;

        for (uint16_t p = 0; p < StrawId::_nupanels; ++p)
            for (size_t l_idx = 0; l_idx < panel_dof_labels[p].size(); ++l_idx)
                label_info_file << "Panel," << p << "," << l_idx << "," << panel_dof_labels[p][l_idx]
                    << std::endl;

        label_info_file.close();

        std::cout << "AlignTrackCollector: wrote labels to " << _labels_filename << std::endl;
    }

    void AlignTrackCollector::beginRun(art::Run const&) { return; }

    void AlignTrackCollector::endJob()
    {
        // ensure the file is closed once the job finishes
        millepede->~Mille();

        if (_diag > 0)
            std::cout << "AlignTrackCollector: wrote " << tracks_written << " tracks to "
                << _output_filename << std::endl;
    }

    bool AlignTrackCollector::filter_CosmicTrackSeedCollection(art::Event const& event,
            Tracker const& tracker,
            StrawResponse const& _srep,
            CosmicTrackSeedCollection const& coscol)
    {
        bool wrote_track = false; // did we write any tracks at all?
        std::cout << "called by event processor" << std::endl;
        // dedicated to CosmicTrackSeedCollection
        for (CosmicTrackSeed const& sts : coscol)
        {
            CosmicTrack const& st = sts._track;
            TrkFitFlag const& status = sts._status;

            std::cout << "found track object" << std::endl;

            if (!status.hasAllProperties(TrkFitFlag::helixOK))
            {
                continue;
            }

            if (_diag > 0 && st.converged)
            {
                // fill seed diagnostics
                for (ComboHit const& hit : sts.hits())
                {
                    plane_seedtracks->Fill(hit.strawId().plane());
                } // Ask Richie about seed chi squared
            }

            if (!st.converged || (!st.minuit_converged && !seed_only))
            {
                continue;
            }


            if (!seed_only && isnan(st.MinuitFitParams.A0))
            {
                continue;
            }

            std::cout << "Have track." << std::endl;

            XYZVec track_pos(st.MinuitFitParams.A0, st.MinuitFitParams.B0, 0);
            XYZVec track_dir(st.MinuitFitParams.A1, st.MinuitFitParams.B1, 1);

            if (seed_only)
            {
                track_pos = st.FitEquationXYZ.Pos;
                track_dir = st.FitEquationXYZ.Dir;
            }
            double A0 = track_pos.X();
            double A1 = track_dir.X();
            double B0 = track_pos.Y();
            double B1 = track_dir.Y();

            double chisq = 0;
            double chsqndof = 0;

            bool wrote_hits = false; // did we write any hits for the track?

            // get residuals and their derivatives with respect
            // to all local and global parameters
            // get also plane id hit by straw hits
            for (ComboHit const& straw_hit : sts._straw_chits)
            {
                // straw and plane info
                StrawId const& straw_id = straw_hit.strawId();
                Straw const& straw = tracker.getStraw(straw_id);
                auto plane_id = straw_id.plane();
                auto panel_uid = straw_id.uniquePanel();
                auto panel_id = straw_id.getPanelId();

                // geometry info
                auto const& plane_origin = tracker.getPlane(plane_id).origin();

                // TODO: check consistency with AlignedTrackerMaker rotation pivot
                auto const& panel_origin = tracker.getPanel(panel_id).straw0MidPoint();
                auto const& straw_mp = straw.getMidPoint();
                auto const& wire_dir = straw.getDirection().unit();

                // now calculate the derivatives.
                auto derivativesLocal = CosmicTrack_DCA_LocalDeriv(
                        A0, B0, A1,
                        B1, straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(),
                        wire_dir.y(), wire_dir.z(),

                        // warning: this is not changed in AlignedTrackerMaker
                        // this is a problem if our starting geometry is not simply the
                        // nominal geometry
                        plane_origin.x(), plane_origin.y(), plane_origin.z(), panel_origin.x(),
                        panel_origin.y(), panel_origin.z());

                auto derivativesGlobal = CosmicTrack_DCA_GlobalDeriv(
                        A0, B0, A1,
                        B1, straw_mp.x(), straw_mp.y(), straw_mp.z(), wire_dir.x(),
                        wire_dir.y(), wire_dir.z(),

                        // warning: this is not changed in AlignedTrackerMaker
                        // this is a problem if our starting geometry is not simply the
                        // nominal geometry
                        plane_origin.x(), plane_origin.y(), plane_origin.z(), panel_origin.x(),
                        panel_origin.y(), panel_origin.z());

                Hep3Vector td(A1, B1, 1);
                td = td.unit();
                Hep3Vector rperp = td - (td.dot(straw.getDirection()) * straw.getDirection());

                // Distance of Closest Approach (DOCA)
                TwoLinePCA_XYZ PCA = TwoLinePCA_XYZ(track_pos, track_dir, Geom::toXYZVec(straw_mp),
                        Geom::toXYZVec(wire_dir), 1.e-8);

                double phi = rperp.theta();
                double drift_distance = // TODO: check this with Dave
                    _srep.driftTimeToDistance(straw_id, straw_hit.driftTime(),
                                             phi);

                // signed DCA.
                double residual = (PCA.LRambig() * PCA.dca()) - drift_distance;
                double residual_error =
                    _srep.driftDistanceError(straw_id, drift_distance, phi, PCA.dca());

                std::cout << "residual " << residual << " +- " << residual_error << std::endl;

                if (isnan(residual))
                    continue;

                chisq += residual * residual;
                chsqndof++;

                // TODO: this is horrendous. surely there is a better way to combine the DOF labels
                std::vector<int> global_dof_labels;

                global_dof_labels.reserve(_dof_per_plane + _dof_per_panel);

                // Convention note:
                // The DoF order is : (planes) dx, dy, dz, a, b, g, followed by (panel) dx, dy, dz, dz,
                // a, b, g This is reflected also in the generated DOCA derivatives.
                global_dof_labels.insert(global_dof_labels.end(), plane_dof_labels[plane_id].begin(),
                        plane_dof_labels[plane_id].end());
                global_dof_labels.insert(global_dof_labels.end(), panel_dof_labels[panel_uid].begin(),
                        panel_dof_labels[panel_uid].end());
                if (_diag > 0 && global_dof_labels.size() != 12) std::cout << "WARNING: we should have 12 labels!" << std::endl;
                // write the hit to the track buffer
                millepede->mille(derivativesLocal.size(), derivativesLocal.data(),
                        derivativesGlobal.size(), derivativesGlobal.data(),
                        global_dof_labels.data(), (float)residual, (float)residual_error);
                wrote_hits = true;
                // diagnostic information
                if (_diag > 0)
                {
                    residuum->Fill(residual);
                    resid_err->Fill(residual_error);
                    plane_tracks->Fill(straw_id.plane());
                    plane_residsum->Fill(straw_id.plane(), residual);

                    drift_time->Fill(straw_hit.driftTime());
                    drift_dist->Fill(drift_distance);
                    doca_h->Fill((PCA.LRambig() * PCA.dca()));

                }
            }
            if (wrote_hits)
            {
                tracks_written++;

                if (chsqndof > 4 && _diag > 0)
                {

                    chsqndof -= 4.0; // 4 track parameters
                    boost::math::chi_squared mydist(chsqndof);
                    track_pvalue->Fill(boost::math::cdf(mydist, chisq));

                    track_chisq->Fill(chisq / chsqndof);
                }
                // Write the track buffer to file if we wrote a track
                millepede->end();

                wrote_track = true;
            }
        }
        return wrote_track;
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

        auto stH = event.getValidHandle<CosmicTrackSeedCollection>(_costag);
        if (stH.product() == 0) return;

        CosmicTrackSeedCollection const& coscol = *stH.product();

        filter_CosmicTrackSeedCollection(event, tracker, _srep, coscol);
    }

}; // namespace mu2e

using mu2e::AlignTrackCollector;
DEFINE_ART_MODULE(AlignTrackCollector);
