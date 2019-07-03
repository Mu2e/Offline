// Andrei Gaponenko, 2016

#include "CutAndCountAnalysis/inc/TrackLevelCuts.hh"

#include <string>
#include <vector>
#include <iostream>

#include "art_root_io/TFileDirectory.h"
#include "art_root_io/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"

#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib_except/exception.h"

namespace mu2e {
  namespace CutAndCount {
    //----------------------------------------------------------------
    // use X-macros to maintain bin labels of cut flow histograms
    // in sync with the cut list
#undef X
#define TRACK_LEVEL_CUTS                        \
    X(status)                                   \
    X(quality)                                  \
    X(pitch)                                    \
    X(d0)                                       \
    X(maxd)                                     \
    X(t0)                                       \
    X(caloMatch)                                \
    X(caloMatchChi2)                            \
    X(caloClusterEnergy)                        \
    X(particleID)                               \
    X(momentum)                                 \
    X(accepted)

    enum class TrkCutNumber {
#define X(entry) entry,
      TRACK_LEVEL_CUTS
      CUTS_END
#undef X
    };

    namespace {
      void set_cut_bin_labels(TAxis* ax) {
#define X(entry) ax->SetBinLabel(1 + int(TrkCutNumber::entry), #entry);
        TRACK_LEVEL_CUTS
#undef X
          }
    }

#undef TRACK_LEVEL_CUTS

  } // namespace CutAndCount

  using namespace CutAndCount;

  art::TFileDirectory TrackLevelCuts::getdir(art::TFileDirectory orig, const std::string& relpath) {
    return relpath.empty() ? orig : orig.mkdir(relpath);
  }

  TrackLevelCuts::TrackLevelCuts(const fhicl::ParameterSet& pset, art::TFileDirectory tfdir, const std::string histdir)
    : TrackLevelCuts(fhicl::Table<PhysicsCuts>{pset, std::set<std::string>{}}(), tfdir, histdir)
  {}

  TrackLevelCuts::TrackLevelCuts(const PhysicsCuts& pc, art::TFileDirectory tfdir, const std::string histdir)
    : TrackLevelCuts(pc, getdir(tfdir,histdir))
  {}

  TrackLevelCuts::TrackLevelCuts(const PhysicsCuts& pc, art::TFileDirectory tf)
    : cuts_(pc)
    , pid_dt_(pc.caloCutsEnabled() ? new PIDdt{pc.caloCuts().pid_dt_conf()} : nullptr)
    , pid_ep_(pc.caloCutsEnabled() ? new PIDEp{pc.caloCuts().pid_ep_conf()} : nullptr)
    , h_cuts_p_{tf.make<TH1D>("cuts_p", "Unweighted events before cut", double(TrkCutNumber::CUTS_END), -0.5, double(TrkCutNumber::CUTS_END)-0.5)}
    , h_cuts_r_{tf.make<TH1D>("cuts_r", "Unweighted events rejected by cut", double(TrkCutNumber::CUTS_END), -0.5, double(TrkCutNumber::CUTS_END)-0.5)}
    , w_cuts_p_{tf.make<TH1D>("wcuts_p", "Weighted events before cut", double(TrkCutNumber::CUTS_END), -0.5, double(TrkCutNumber::CUTS_END)-0.5)}
    , w_cuts_r_{tf.make<TH1D>("wcuts_r", "Weighted events rejected by cut", double(TrkCutNumber::CUTS_END), -0.5, double(TrkCutNumber::CUTS_END)-0.5)}
    , trkqual_{tf.make<TH1D>("trkqual", "trkqual before cut", 100, 0., 2.)}
    , td_{tf.make<TH1D>("td", "Track tan(lambda) before cut", 100, 0.5, 1.5)}
    , d0_{tf.make<TH1D>("d0", "Track d0 before cut", 300, -150., +150.)}
    , rmax_{tf.make<TH1D>("rmax", "Track d0+2/om  before cut", 120, 300., 900.)}
    , t0_{tf.make<TH1D>("t0", "Track t0  before cut", 170, 0., 1700.)}
    , caloMatchChi2_{tf.make<TH1D>("caloMatchCHi2", "Calo match chi2 before cut", 100, 0., 300.)}
    , caloClusterEnergy_{tf.make<TH1D>("caloClusterEnergy", "Calo cluster energy before cut", 150, 0., 150.)}
    , pidVariable_{tf.make<TH1D>("particleID", "PID variable before cut", 1000, -50., 50.)}
    , momentum_{tf.make<TH1D>("momentum", "Track momentum  before cut", 1000, 98., 108.)}
  {
    using CutAndCount::TrkCutNumber;
    using CutAndCount::set_cut_bin_labels;

    h_cuts_p_->SetStats(kFALSE);
    set_cut_bin_labels(h_cuts_p_->GetXaxis());
    h_cuts_p_->SetOption("hist text");

    h_cuts_r_->SetStats(kFALSE);
    set_cut_bin_labels(h_cuts_r_->GetXaxis());
    h_cuts_r_->SetOption("hist text");

    w_cuts_p_->SetStats(kFALSE);
    set_cut_bin_labels(w_cuts_p_->GetXaxis());
    w_cuts_p_->SetOption("hist text");
    w_cuts_p_->Sumw2();

    w_cuts_r_->SetStats(kFALSE);
    set_cut_bin_labels(w_cuts_r_->GetXaxis());
    w_cuts_r_->SetOption("hist text");
    w_cuts_r_->Sumw2();

    trkqual_->Sumw2();
    td_->Sumw2();
    d0_->Sumw2();
    rmax_->Sumw2();
    t0_->Sumw2();
    caloMatchChi2_->Sumw2();
    caloClusterEnergy_->Sumw2();
    pidVariable_->Sumw2();
    momentum_->Sumw2();
  }

  //================================================================
  bool TrackLevelCuts::accepted(const art::Ptr<KalRep>& trk, const art::Event& event, const KalDiag& kdiag, const EventWeightHelper& wh) {
    using CutAndCount::TrkCutNumber;

    TrkCutNumber c = processTrack(trk, event, kdiag, wh);
    h_cuts_r_->Fill(double(c));
    w_cuts_r_->Fill(double(c), wh.weight());

    for(int cut=0; cut<=int(c); cut++) {
      h_cuts_p_->Fill(cut);
      w_cuts_p_->Fill(cut, wh.weight());
    }

    return c==TrkCutNumber::accepted;
  }

  //================================================================
  CutAndCount::TrkCutNumber TrackLevelCuts::processTrack(const art::Ptr<KalRep>& trk,
                                                         const art::Event& evt,
                                                         const KalDiag& kdiag,
                                                         const EventWeightHelper& wh)
  {
    using CutAndCount::TrkCutNumber;

    if(!trk->fitCurrent()) {
      throw cet::exception("BADINPUT")<<"TrackLevelCuts: do not know what to do with a fitCurrent==0 track\n";
    }

    TrkInfo track;
    kdiag.fillTrkInfo(trk.get(), track);
    TrkFitInfo track_ent;
    kdiag.fillTrkFitInfo(trk.get(), track_ent);

    // _status > 0 the means code found a track
    if(track._status <= 0) {
      return TrkCutNumber::status;
    }

    trkqual_->Fill(track._trkqual, wh.weight());
    if((track._trkqual < cuts_.trkqual()) ||
       ((cuts_.trkqualmax() > 0.) && (cuts_.trkqualmax() < track._trkqual) ) ) {
      return TrkCutNumber::quality;
    }

    const helixpar& th = track_ent._fitpar;
    td_->Fill(th._td, wh.weight());
    if((th._td < cuts_.tdmin())||(th._td > cuts_.tdmax())) {
      return TrkCutNumber::pitch;
    }

    d0_->Fill(th._d0, wh.weight());
    if((th._d0 < cuts_.d0min())||(th._d0 > cuts_.d0max())) {
      return TrkCutNumber::d0;
    }

    const double maxd = th._d0 + 2./th._om;
    rmax_->Fill(maxd, wh.weight());
    if((maxd < cuts_.mdmin())||(maxd > cuts_.mdmax())) {
      return TrkCutNumber::maxd;
    }

    t0_->Fill(track._t0, wh.weight());
    if((track._t0 < cuts_.t0min())||(track._t0 > cuts_.t0max())) {
      return TrkCutNumber::t0;
    }

    //----------------------------------------------------------------
    // Here we start using calorimeter info

    if(cuts_.caloCutsEnabled()) {
      const TrackCaloCuts& ccuts = cuts_.caloCuts();

      const auto* cm = findCaloMatch(trk, evt);
      if(!cm) {
        return TrkCutNumber::caloMatch;
      }

      caloMatchChi2_->Fill(cm->chi2(), wh.weight());
      if(cm->chi2() > ccuts.matchChi2()) {
        return TrkCutNumber::caloMatchChi2;
      }

      const double clusterEnergy = cm->caloCluster()->energyDep();
      caloClusterEnergy_->Fill(clusterEnergy, wh.weight());
      if((clusterEnergy < ccuts.emin())||(ccuts.emax() < clusterEnergy)) {
        return TrkCutNumber::caloClusterEnergy;
      }

      double pidvar = pid_dt_->value(cm->dt());
      // value at cutoff() means the inputs are not consistent with signal
      if(pidvar > pid_dt_->cutoff()) {
        const double pidvar_ep = pid_ep_->value(cm->ep(), cm->ds());
        if(pidvar_ep > pid_ep_->cutoff()) {
          pidvar += pidvar_ep;
        }
        else {
          pidvar = pid_ep_->cutoff();
        }
      }

      pidVariable_->Fill(pidvar, wh.weight());
      if(pidvar < ccuts.pidCut()) {
        return TrkCutNumber::particleID;
      }

    }

    //----------------------------------------------------------------
    // The analysis momentum window cut
    const double fitmom = track_ent._fitmom;
    momentum_->Fill(fitmom, wh.weight());
    if((fitmom < cuts_.pmin())||(fitmom > cuts_.pmax())) {
      return TrkCutNumber::momentum;
    }

    return TrkCutNumber::accepted;

  } // processTrack()

  //================================================================
  const TrackClusterMatch* TrackLevelCuts::findCaloMatch(const art::Ptr<KalRep>& trk, const art::Event& evt) {
    //auto ihmatch = evt.getValidHandle<TrkCaloMatchCollection>(caloMatchDemInput_);
    auto ihmatch = evt.getValidHandle<TrackClusterMatchCollection>(cuts_.caloCuts().clusterInput());
    for(const auto& match: *ihmatch) {
      if(match.textrapol()->trk() == trk) {
        return &match;
      }
    }
    return nullptr;
  }

  //================================================================

} // namespace mu2e
