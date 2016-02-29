// A module to produce a spectrum of reconstructed electrons around
// the signal window, with all the analysis cuts.  This is the common
// part of different background analyses.
//
// Analysis-specific variations are introduced at the job
// configuration level (fcl).  For example, RPC analysis needs to
// weight events according to the proper time of the stopped pion.
// The weight should be computed in a separate module and used as in
// input here.
//
// This module
//
//    - is not supposed to have any MC truth dependencies.  There is
//      an optional EventWeight input that may be computed from MC
//      truth, but weight computation is outside of this module.
//
//    - is not meant for detailed studies.  Studies to come up with
//      cut definitions should be done elsewhere, not added to this
//      code.  Here we just plug a set of pre-define "standard" cuts
//      and get an answer (track count in the signal region).
//
// Andrei Gaponenko, 2016

#include <string>
#include <vector>
#include <iostream>

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Utilities/InputTag.h"

#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"


#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib/exception.h"

#include "boost/noncopyable.hpp"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "TrkDiag/inc/KalDiag.hh"

#include "TrackCaloMatching/inc/TrackClusterMatch.hh"
//#include "RecoDataProducts/inc/TrkCaloMatchCollection.hh"
//#include "RecoDataProducts/inc/TrkCaloMatch.hh"
#include "RecoDataProducts/inc/TrkCaloIntersect.hh"

#include "TH1.h"
#include "TH2.h"

namespace mu2e {
  namespace {
    //----------------------------------------------------------------
    // Histograms filled during track selection: neither "all" nor "accepted" tracks
    struct TrkCutHist: private boost::noncopyable {
      explicit TrkCutHist(art::TFileDirectory tfdir, const std::string& relpath);
      explicit TrkCutHist(art::TFileDirectory tfdir);

      TH1 *trkqual;
      TH1 *td;
      TH1 *d0;
      TH1 *rmax;
      TH1 *t0;
      TH1 *momentum;
    private:
      art::TFileDirectory getdir(art::TFileDirectory orig, const std::string& relpath);
    };

    art::TFileDirectory TrkCutHist::getdir(art::TFileDirectory orig, const std::string& relpath) {
      return relpath.empty() ? orig : orig.mkdir(relpath);
    }

    TrkCutHist::TrkCutHist(art::TFileDirectory tfdir, const std::string& relpath)
      : TrkCutHist(getdir(tfdir, relpath))
    {}

    TrkCutHist::TrkCutHist(art::TFileDirectory tf)
      : trkqual{tf.make<TH1D>("trkqual", "trkqual before cut", 100, 0., 1.)}
      , td{tf.make<TH1D>("td", "Track tan(lambda) before cut", 100, 0.5, 1.5)}
      , d0{tf.make<TH1D>("d0", "Track d0 before cut", 300, -150., +150.)}
      , rmax{tf.make<TH1D>("rmax", "Track d0+2/om  before cut", 120, 300., 900.)}
      , t0{tf.make<TH1D>("t0", "Track t0  before cut", 170, 0., 1700.)}
      , momentum{tf.make<TH1D>("momentum", "Track momentum  before cut", 500, 98., 108.)}
    {}

  } // namespace{}

  class CutAndCountAnalysis : public art::EDAnalyzer {
  public:

    explicit CutAndCountAnalysis(const fhicl::ParameterSet& pset);
    virtual void analyze(const art::Event& event) override;

  private:
    // A structure to hold values of physics cuts
    struct PhysicsCuts {
      double trkqual;
      // tan(lambda)
      double tdmin;
      double tdmax;
      // closest aproach to (0,0)
      double d0min;
      double d0max;
      // max distance from (0,0)
      double mdmin;
      double mdmax;
      // time
      double t0min;
      double t0max;
      // momentum window
      double pmin;
      double pmax;

      PhysicsCuts(const fhicl::ParameterSet& pset)
        : trkqual(pset.get<double>("trkqual"))
        , tdmin(pset.get<double>("tdmin"))
        , tdmax(pset.get<double>("tdmax"))
        , d0min(pset.get<double>("d0min"))
        , d0max(pset.get<double>("d0max"))
        , mdmin(pset.get<double>("mdmin"))
        , mdmax(pset.get<double>("mdmax"))
        , t0min(pset.get<double>("t0min"))
        , t0max(pset.get<double>("t0max"))
        , pmin(pset.get<double>("pmin"))
        , pmax(pset.get<double>("pmax"))
      {}
    };

    //----------------------------------------------------------------
    typedef std::vector<art::InputTag> InputTags;
    art::InputTag trackDemInput_;
    art::InputTag caloMatchDemInput_;
    InputTags weights_;
    PhysicsCuts cuts_;
    KalDiag kdiag_;

    TH1 *h_cuts_p_;
    TH1 *h_cuts_r_;
    TH1 *hNumAcceptedTracks_;
    TrkCutHist hTrkCuts_;

    //----------------------------------------------------------------
    // use X-macros to maintain bin labels of cut flow histograms
    // in sync with the cut list
#undef X
#define TRACK_LEVEL_CUTS                              \
    X(status)                                         \
    X(quality)                                        \
    X(pitch)                                          \
    X(d0)                                             \
    X(maxd)                                           \
    X(t0)                                             \
    X(caloMatch)                                      \
    X(momentum)                                       \
    X(accepted)

    enum class TrkCut {
#define X(entry) entry,
      TRACK_LEVEL_CUTS
      CUTS_END
#undef X
    };

    void set_cut_bin_labels(TAxis* ax) {
#define X(entry) ax->SetBinLabel(1 + int(TrkCut::entry), #entry);
      TRACK_LEVEL_CUTS
#undef X
    }

#undef TRACK_LEVEL_CUTS

    TrkCut processTrack(const art::Ptr<KalRep>& trk, const art::Event& evt);
    const TrackClusterMatch* findCaloMatch(const art::Ptr<KalRep>& trk, const art::Event& evt);
  };

  //================================================================
  CutAndCountAnalysis::CutAndCountAnalysis(const fhicl::ParameterSet& pset)
    : art::EDAnalyzer(pset)
    , trackDemInput_(pset.get<art::InputTag>("trackDemInput"))
    , caloMatchDemInput_(pset.get<art::InputTag>("caloMatchDemInput"))
    , weights_(pset.get<std::vector<art::InputTag> >("weights"))
    , cuts_(pset.get<fhicl::ParameterSet>("physicsCuts"))
    , kdiag_(pset.get<fhicl::ParameterSet>("kalDiag"))
    , hTrkCuts_(*art::ServiceHandle<art::TFileService>(), "trkcuts")
  {
    art::ServiceHandle<art::TFileService> tfs;
    h_cuts_p_ = tfs->make<TH1D>("cuts_p", "Events before cut", double(TrkCut::CUTS_END), -0.5, double(TrkCut::CUTS_END)-0.5);
    h_cuts_p_->SetStats(kFALSE);
    set_cut_bin_labels(h_cuts_p_->GetXaxis());
    h_cuts_p_->SetOption("hist text");

    h_cuts_r_ = tfs->make<TH1D>("cuts_r", "Events rejected by cut", double(TrkCut::CUTS_END), -0.5, double(TrkCut::CUTS_END)-0.5);
    h_cuts_r_->SetStats(kFALSE);
    set_cut_bin_labels(h_cuts_r_->GetXaxis());
    h_cuts_r_->SetOption("hist text");

    hNumAcceptedTracks_ = tfs->make<TH1D>("numAcceptedTracks", "Number of accepted tracks per event", 5, -0.5, 4.5);
  }

  //================================================================
  void CutAndCountAnalysis::analyze(const art::Event& event) {

    auto ih = event.getValidHandle<KalRepPtrCollection>(trackDemInput_);

    int acceptedTracksCount = 0;
    for(const auto& ptr: *ih) {
      TrkCut c = processTrack(ptr, event);
      h_cuts_r_->Fill(double(c));
      for(int cut=0; cut<=int(c); cut++) {
        h_cuts_p_->Fill(cut);
      }
      if(c==TrkCut::accepted) {
        ++acceptedTracksCount;
      }
    }

    hNumAcceptedTracks_->Fill(acceptedTracksCount);
  }

  //================================================================
  CutAndCountAnalysis::TrkCut CutAndCountAnalysis::processTrack(const art::Ptr<KalRep>& trk, const art::Event& evt) {
    if(!trk->fitCurrent()) {
      throw cet::exception("BADINPUT")<<"CutAndCountAnalysis: do not know what to do with a fitCurrent==0 track\n";
    }

    TrkInfo track;
    kdiag_.fillTrkInfo(trk.get(), track);

    if(track._fitstatus != 1) {
      return TrkCut::status;
    }

    hTrkCuts_.trkqual->Fill(track._trkqual);
    if(track._trkqual < cuts_.trkqual) {
      return TrkCut::quality;
    }

    const helixpar& th = track._ent._fitpar;
    hTrkCuts_.td->Fill(th._td);
    if((th._td < cuts_.tdmin)||(th._td > cuts_.tdmax)) {
      return TrkCut::pitch;
    }

    hTrkCuts_.d0->Fill(th._d0);
    if((th._d0 < cuts_.d0min)||(th._d0 > cuts_.d0max)) {
      return TrkCut::d0;
    }

    const double maxd = th._d0 + 2./th._om;
    hTrkCuts_.rmax->Fill(maxd);
    if((maxd < cuts_.mdmin)||(maxd > cuts_.mdmax)) {
      return TrkCut::maxd;
    }

    hTrkCuts_.t0->Fill(track._t0);
    if((track._t0 < cuts_.t0min)||(track._t0 > cuts_.t0max)) {
      return TrkCut::t0;
    }

    //----------------------------------------------------------------
    // Here we start using calorimeter info
    const auto* cm = findCaloMatch(trk, evt);
    if(!cm) {
      return TrkCut::caloMatch;
    }

    //----------------------------------------------------------------
    // The analysis momentum window cut
    const double fitmom = track._ent._fitmom;
    hTrkCuts_.momentum->Fill(fitmom);
    if((fitmom < cuts_.pmin)||(fitmom > cuts_.pmax)) {
      return TrkCut::momentum;
    }

    return TrkCut::accepted;

  } // processTrack()

  //================================================================
  const TrackClusterMatch* CutAndCountAnalysis::findCaloMatch(const art::Ptr<KalRep>& trk, const art::Event& evt) {
    //auto ihmatch = evt.getValidHandle<TrkCaloMatchCollection>(caloMatchDemInput_);
    auto ihmatch = evt.getValidHandle<TrackClusterMatchCollection>(caloMatchDemInput_);
    std::cout<<"----------------------------------------------------------------"<<std::endl;
    for(const auto& match: *ihmatch) {
      std::cout<<"Looking at match.trk = "<<match.textrapol()->trk()
               <<", our trk = "<<trk
               <<std::endl;
      if(match.textrapol()->trk() == trk) {
        return &match;
      }
    }

    return nullptr;
  }

  //================================================================

} // namespace mu2e

DEFINE_ART_MODULE(mu2e::CutAndCountAnalysis);
