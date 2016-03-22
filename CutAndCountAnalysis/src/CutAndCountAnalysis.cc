// Andrei Gaponenko, 2016

#include "CutAndCountAnalysis/inc/CutAndCountAnalysis.hh"

#include <string>
#include <vector>
#include <iostream>

#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Registry/ServiceHandle.h"


#include "CLHEP/Units/SystemOfUnits.h"

#include "cetlib/exception.h"

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
    X(momentum)                                 \
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

  } // namespace CutAndCount

  CutAndCountAnalysis::CutAndCountAnalysis(const fhicl::ParameterSet& pset, art::TFileDirectory tfdir)
    : trackDemInput_(pset.get<art::InputTag>("trackDemInput"))
    , caloMatchDemInput_(pset.get<art::InputTag>("caloMatchDemInput"))
    , cuts_(pset.get<fhicl::ParameterSet>("physicsCuts"))
    , wh_(pset.get<fhicl::ParameterSet>("weight"), *art::ServiceHandle<art::TFileService>(), "weight")
    , kdiag_(pset.get<fhicl::ParameterSet>("kalDiag"))
    , hTrkCuts_(tfdir, "trkcuts")
  {
    using CutAndCount::TrkCut;
    using CutAndCount::set_cut_bin_labels;

    h_cuts_p_ = tfdir.make<TH1D>("cuts_p", "Unweighted events before cut", double(TrkCut::CUTS_END), -0.5, double(TrkCut::CUTS_END)-0.5);
    h_cuts_p_->SetStats(kFALSE);
    set_cut_bin_labels(h_cuts_p_->GetXaxis());
    h_cuts_p_->SetOption("hist text");

    h_cuts_r_ = tfdir.make<TH1D>("cuts_r", "Unweighted events rejected by cut", double(TrkCut::CUTS_END), -0.5, double(TrkCut::CUTS_END)-0.5);
    h_cuts_r_->SetStats(kFALSE);
    set_cut_bin_labels(h_cuts_r_->GetXaxis());
    h_cuts_r_->SetOption("hist text");

    w_cuts_p_ = tfdir.make<TH1D>("wcuts_p", "Weighted events before cut", double(TrkCut::CUTS_END), -0.5, double(TrkCut::CUTS_END)-0.5);
    w_cuts_p_->SetStats(kFALSE);
    set_cut_bin_labels(w_cuts_p_->GetXaxis());
    w_cuts_p_->SetOption("hist text");
    w_cuts_p_->Sumw2();

    w_cuts_r_ = tfdir.make<TH1D>("wcuts_r", "Weighted events rejected by cut", double(TrkCut::CUTS_END), -0.5, double(TrkCut::CUTS_END)-0.5);
    w_cuts_r_->SetStats(kFALSE);
    set_cut_bin_labels(w_cuts_r_->GetXaxis());
    w_cuts_r_->SetOption("hist text");
    w_cuts_r_->Sumw2();

    hNumAcceptedTracks_ = tfdir.make<TH1D>("numAcceptedTracks", "Number of accepted tracks per event (weighted)", 5, -0.5, 4.5);
    hNumAcceptedTracks_->Sumw2();
  }

  //================================================================
  bool CutAndCountAnalysis::accepted(const art::Event& event) {
    using CutAndCount::TrkCut;

    bool passed = false;

    wh_.update(event);

    auto ih = event.getValidHandle<KalRepPtrCollection>(trackDemInput_);

    int acceptedTracksCount = 0;
    for(const auto& ptr: *ih) {
      TrkCut c = processTrack(ptr, event);
      h_cuts_r_->Fill(double(c));
      w_cuts_r_->Fill(double(c), wh_.weight());

      for(int cut=0; cut<=int(c); cut++) {
        h_cuts_p_->Fill(cut);
        w_cuts_p_->Fill(cut, wh_.weight());
      }
      if(c==TrkCut::accepted) {
        ++acceptedTracksCount;
        passed = true;
      }
    }

    hNumAcceptedTracks_->Fill(acceptedTracksCount, wh_.weight());
    return passed;
  }

  //================================================================
  CutAndCount::TrkCut CutAndCountAnalysis::processTrack(const art::Ptr<KalRep>& trk, const art::Event& evt) {
    using CutAndCount::TrkCut;

    if(!trk->fitCurrent()) {
      throw cet::exception("BADINPUT")<<"CutAndCountAnalysis: do not know what to do with a fitCurrent==0 track\n";
    }

    TrkInfo track;
    kdiag_.fillTrkInfo(trk.get(), track);

    if(track._fitstatus != 1) {
      return TrkCut::status;
    }

    hTrkCuts_.trkqual->Fill(track._trkqual, wh_.weight());
    if(track._trkqual < cuts_.trkqual) {
      return TrkCut::quality;
    }

    const helixpar& th = track._ent._fitpar;
    hTrkCuts_.td->Fill(th._td, wh_.weight());
    if((th._td < cuts_.tdmin)||(th._td > cuts_.tdmax)) {
      return TrkCut::pitch;
    }

    hTrkCuts_.d0->Fill(th._d0, wh_.weight());
    if((th._d0 < cuts_.d0min)||(th._d0 > cuts_.d0max)) {
      return TrkCut::d0;
    }

    const double maxd = th._d0 + 2./th._om;
    hTrkCuts_.rmax->Fill(maxd, wh_.weight());
    if((maxd < cuts_.mdmin)||(maxd > cuts_.mdmax)) {
      return TrkCut::maxd;
    }

    hTrkCuts_.t0->Fill(track._t0, wh_.weight());
    if((track._t0 < cuts_.t0min)||(track._t0 > cuts_.t0max)) {
      return TrkCut::t0;
    }

    //----------------------------------------------------------------
    // Here we start using calorimeter info
    const auto* cm = findCaloMatch(trk, evt);
    if(!cm) {
      return TrkCut::caloMatch;
    }

    hTrkCuts_.caloMatchChi2->Fill(cm->chi2(), wh_.weight());
    if(cm->chi2() > cuts_.caloMatchChi2) {
      return TrkCut::caloMatchChi2;
    }

    const double clusterEnergy = cm->caloCluster()->energyDep();
    hTrkCuts_.caloClusterEnergy->Fill(clusterEnergy, wh_.weight());
    if((clusterEnergy < cuts_.caloemin)||(cuts_.caloemax < clusterEnergy)) {
      return TrkCut::caloClusterEnergy;
    }

    //----------------------------------------------------------------
    // The analysis momentum window cut
    const double fitmom = track._ent._fitmom;
    hTrkCuts_.momentum->Fill(fitmom, wh_.weight());
    if((fitmom < cuts_.pmin)||(fitmom > cuts_.pmax)) {
      return TrkCut::momentum;
    }

    return TrkCut::accepted;

  } // processTrack()

  //================================================================
  const TrackClusterMatch* CutAndCountAnalysis::findCaloMatch(const art::Ptr<KalRep>& trk, const art::Event& evt) {
    //auto ihmatch = evt.getValidHandle<TrkCaloMatchCollection>(caloMatchDemInput_);
    auto ihmatch = evt.getValidHandle<TrackClusterMatchCollection>(caloMatchDemInput_);
    for(const auto& match: *ihmatch) {
      if(match.textrapol()->trk() == trk) {
        return &match;
      }
    }
    return nullptr;
  }

  //================================================================

} // namespace mu2e
