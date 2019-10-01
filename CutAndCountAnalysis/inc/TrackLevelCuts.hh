// Code to produce a spectrum of reconstructed electrons around
// the signal window, with all the analysis cuts.  This is the common
// part of different background analyses.
//
// Analysis-specific variations are introduced at the job
// configuration level (fcl).  For example, RPC analysis needs to
// weight events according to the proper time of the stopped pion.
// The weight should be computed in a separate module and used as in
// input here.
//
// This code
//
//    - is not supposed to have any MC truth dependencies.  There is
//      an optional EventWeight input that may be computed from MC
//      truth, but weight computation is outside of this module.
//
//    - is not meant for detailed studies.  Studies to come up with
//      cut definitions should be done elsewhere, not added to this
//      code.  Here we just plug a set of pre-defined "standard" cuts
//      and get an answer (track count in the signal region).
//
// Andrei Gaponenko, 2016

#ifndef CutAndCountAnalysis_inc_TrackLevelCuts_hh
#define CutAndCountAnalysis_inc_TrackLevelCuts_hh

#include <string>
#include <vector>
#include <memory>

#include "boost/noncopyable.hpp"

#include "art_root_io/TFileDirectory.h"

#include "fhiclcpp/types/Atom.h"
#include "fhiclcpp/types/Sequence.h"
#include "fhiclcpp/types/Table.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "canvas/Utilities/InputTag.h"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "TrkDiag/inc/KalDiag.hh"

#include "Mu2eUtilities/inc/EventWeightHelper.hh"

#include "RecoDataProducts/inc/TrackClusterMatch.hh"
//#include "RecoDataProducts/inc/TrkCaloMatchCollection.hh"
//#include "RecoDataProducts/inc/TrkCaloMatch.hh"
#include "RecoDataProducts/inc/TrkCaloIntersect.hh"

#include "ParticleID/inc/PIDLogLRatio.hh"
#include "ParticleID/inc/PIDLogL1D.hh"
#include "ParticleID/inc/PIDLogLEp.hh"

#include "TH1.h"
#include "TH2.h"

namespace mu2e {

  namespace CutAndCount {
    enum class TrkCutNumber;
  }

  class TrackLevelCuts: private boost::noncopyable {
  public:
    typedef PIDLogLRatio<PIDLogL1D> PIDdt;
    typedef PIDLogLRatio<PIDLogLEp> PIDEp;

    // A structure to hold a subset of physics cuts related to track-calo matching
    struct TrackCaloCuts {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<art::InputTag> clusterInput{ Name("clusterInput"), Comment("Tag of the calo cluster collection to match to the tracks.")};

      fhicl::Atom<double> matchChi2{ Name("matchChi2"), Comment("High cut on chi2 of track-calo match") };

      fhicl::Atom<double> emin{Name("emin"), Comment("Min energy of matched calo cluster")};

      fhicl::Atom<double> emax{Name("emax"), Comment("Max energy of matched calo cluster")};

      fhicl::Table<PIDdt::Config> pid_dt_conf{Name("PIDdt"), Comment("dt based PID config")};
      fhicl::Table<PIDEp::Config> pid_ep_conf{Name("PIDEp"), Comment("E/p based PID config")};

      fhicl::Atom<double> pidCut{Name("pidCut"), Comment("Low cut on total PID variable")};
    };

    // A top level structure to hold values of physics cuts
    struct PhysicsCuts {
      using Name=fhicl::Name;
      using Comment=fhicl::Comment;

      fhicl::Atom<double> trkqual{Name("trkqual"), Comment("Low cut on trkqual") };
      fhicl::Atom<double> trkqualmax{Name("trkqualmax"), Comment("Optional high cut on trkqual. A non-positive value disables the cut."), -1. };

      fhicl::Atom<double> tdmin{Name("tdmin"), Comment("Low cut on tan(dip)")};
      fhicl::Atom<double> tdmax{Name("tdmax"), Comment("High cut on tan(dip)")};
      fhicl::Atom<double> d0min{Name("d0min"), Comment("Low cut on closest track approach to detector axis")};
      fhicl::Atom<double> d0max{Name("d0max"), Comment("High cut on closest track approach to detector axis")};
      fhicl::Atom<double> mdmin{Name("mdmin"), Comment("Low cut on farthest track distance from detector axis")};
      fhicl::Atom<double> mdmax{Name("mdmax"), Comment("High cut on farthest track distance from detector axis")};
      fhicl::Atom<double> t0min{Name("t0min"), Comment("Beginning of the live time gate")};
      fhicl::Atom<double> t0max{Name("t0max"), Comment("End of the live time gate")};

      fhicl::Atom<bool> caloCutsEnabled{ Name("caloCutsEnabled"),
          Comment("Set this to false to turn off calorimeter matching and related cuts. "),
          true
          };

      fhicl::Table<TrackCaloCuts> caloCuts{ Name("caloCuts"),
          Comment("Config for calo matching and related cuts"),
          fhicl::MaybeUseFunction([this](){ return caloCutsEnabled(); })
      };

      fhicl::Atom<double> pmin{Name("pmin"), Comment("Low cut on signal track momentum")};
      fhicl::Atom<double> pmax{Name("pmax"), Comment("High cut on signal track momentum")};
    };

    explicit TrackLevelCuts(const PhysicsCuts& pc, art::TFileDirectory tfdir, const std::string histdir);
    explicit TrackLevelCuts(const fhicl::ParameterSet& pset, art::TFileDirectory tfdir, const std::string histdir);

    bool accepted(const art::Ptr<KalRep>& trk,
                  const art::Event& event,
                  const KalDiag& kdiag,
                  const EventWeightHelper& wh);

  private:
    art::TFileDirectory getdir(art::TFileDirectory orig, const std::string& relpath);
    explicit TrackLevelCuts(const PhysicsCuts& pset, art::TFileDirectory finaldir);

    //----------------------------------------------------------------
    PhysicsCuts cuts_;
    std::unique_ptr<PIDdt> pid_dt_;
    std::unique_ptr<PIDEp> pid_ep_;

    TH1 *h_cuts_p_;
    TH1 *h_cuts_r_;
    TH1 *w_cuts_p_;
    TH1 *w_cuts_r_;

    TH1 *trkqual_;
    TH1 *td_;
    TH1 *d0_;
    TH1 *rmax_;
    TH1 *t0_;
    TH1 *caloMatchChi2_;
    TH1 *caloClusterEnergy_;
    TH1 *pidVariable_;
    TH1 *momentum_;

    CutAndCount::TrkCutNumber processTrack(const art::Ptr<KalRep>& trk, const art::Event& evt, const KalDiag& kdiag, const EventWeightHelper& wh);
    const TrackClusterMatch* findCaloMatch(const art::Ptr<KalRep>& trk, const art::Event& evt);
  };

} // namespace mu2e

#endif/*CutAndCountAnalysis_inc_TrackLevelCuts_hh*/
