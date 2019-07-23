// Andrei Gaponenko, 2016

#include "CutAndCountAnalysis/inc/CutAndCountAnalysis.hh"

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
#define EVENT_LEVEL_CUTS                        \
    X(NoSignalCandidate)                        \
    X(TrackVeto)                                \
    X(accepted)

    enum class EventCutNumber {
#define X(entry) entry,
      EVENT_LEVEL_CUTS
      CUTS_END
#undef X
    };

    namespace {
      void set_cut_bin_labels(TAxis* ax) {
#define X(entry) ax->SetBinLabel(1 + int(EventCutNumber::entry), #entry);
        EVENT_LEVEL_CUTS
#undef X
          }
    }

#undef TRACK_LEVEL_CUTS
  } // namespace CutAndCount

  CutAndCountAnalysis::CutAndCountAnalysis(const fhicl::ParameterSet& pset, art::TFileDirectory tfdir)
    : signalCandidateInput_(pset.get<art::InputTag>("signalTrackInput"))
    , signalTrackCuts_(pset.get<fhicl::ParameterSet>("signalTrackCuts"), tfdir, "signalTrackCuts")
    , wh_(pset.get<fhicl::ParameterSet>("weight"), *art::ServiceHandle<art::TFileService>(), "weight")
    , kdiag_(pset.get<fhicl::ParameterSet>("kalDiag"))
  {
    typedef std::tuple<std::string, art::InputTag, fhicl::ParameterSet> VetoPars;
    auto trackVetoes = pset.get<std::vector<VetoPars> >("trackVetoes");
    vetoDefs_.reserve(trackVetoes.size());
    for(const auto& v: trackVetoes) {
      vetoDefs_.emplace_back(std::get<0>(v),
                             std::get<1>(v),
                             std::make_unique<TrackLevelCuts>(std::get<2>(v), tfdir, std::get<0>(v))
                             );
    }

    //----------------------------------------------------------------
    using CutAndCount::EventCutNumber;
    using CutAndCount::set_cut_bin_labels;

    h_cuts_p_ = tfdir.make<TH1D>("cuts_p", "Unweighted events before cut", double(EventCutNumber::CUTS_END), -0.5, double(EventCutNumber::CUTS_END)-0.5);
    h_cuts_p_->SetStats(kFALSE);
    set_cut_bin_labels(h_cuts_p_->GetXaxis());
    h_cuts_p_->SetOption("hist text");

    h_cuts_r_ = tfdir.make<TH1D>("cuts_r", "Unweighted events rejected by cut", double(EventCutNumber::CUTS_END), -0.5, double(EventCutNumber::CUTS_END)-0.5);
    h_cuts_r_->SetStats(kFALSE);
    set_cut_bin_labels(h_cuts_r_->GetXaxis());
    h_cuts_r_->SetOption("hist text");

    w_cuts_p_ = tfdir.make<TH1D>("wcuts_p", "Weighted events before cut", double(EventCutNumber::CUTS_END), -0.5, double(EventCutNumber::CUTS_END)-0.5);
    w_cuts_p_->SetStats(kFALSE);
    set_cut_bin_labels(w_cuts_p_->GetXaxis());
    w_cuts_p_->SetOption("hist text");
    w_cuts_p_->Sumw2();

    w_cuts_r_ = tfdir.make<TH1D>("wcuts_r", "Weighted events rejected by cut", double(EventCutNumber::CUTS_END), -0.5, double(EventCutNumber::CUTS_END)-0.5);
    w_cuts_r_->SetStats(kFALSE);
    set_cut_bin_labels(w_cuts_r_->GetXaxis());
    w_cuts_r_->SetOption("hist text");
    w_cuts_r_->Sumw2();

    hNumSignalCandidates_ = tfdir.make<TH1D>("numSignalCandidates", "Number of signal candidate tracks per event (weighted) before the vetoes", 5, -0.5, 4.5);
    hNumSignalCandidates_->Sumw2();

    if(!vetoDefs_.empty()) {
      hNumVetoCandidates_ = tfdir.make<TH2D>("numVetoCandidates", "Number of veto tracks per event (weighted) for events with signal candidates", vetoDefs_.size(), -0.5, vetoDefs_.size()-0.5, 5, -0.5, 4.5);
      hNumVetoCandidates_->Sumw2();
      hNumVetoCandidates_->SetOption("text");

      for(unsigned i=0; i<vetoDefs_.size(); ++i) {
        std::string name = std::get<0>(vetoDefs_[i]); // bind to a variable to guarantee lifetime
        hNumVetoCandidates_->GetXaxis()->SetBinLabel(1+i, name.c_str());
      }
    }
  }

  //================================================================
  bool CutAndCountAnalysis::accepted(const art::Event& event) {
    using CutAndCount::EventCutNumber;
    wh_.update(event);

    EventCutNumber c = processEvent(event);
    h_cuts_r_->Fill(double(c));
    w_cuts_r_->Fill(double(c), wh_.weight());

    for(int cut=0; cut<=int(c); cut++) {
      h_cuts_p_->Fill(cut);
      w_cuts_p_->Fill(cut, wh_.weight());
    }

    return (c==EventCutNumber::accepted);
  }

  //================================================================
  CutAndCount::EventCutNumber CutAndCountAnalysis::processEvent(const art::Event& evt) {
    using CutAndCount::EventCutNumber;

    auto ih = evt.getValidHandle<KalRepPtrCollection>(signalCandidateInput_);

    int numSignalCandidates = 0;
    for(const auto& ptr: *ih) {
      if(signalTrackCuts_.accepted(ptr, evt, kdiag_, wh_)) {
        ++numSignalCandidates;
      }
    }

    hNumSignalCandidates_->Fill(numSignalCandidates, wh_.weight());
    if(numSignalCandidates==0) {
      return EventCutNumber::NoSignalCandidate;
    }

    // Impose track vetoes
    for(unsigned iveto=0; iveto < vetoDefs_.size(); ++iveto) {
      auto uemh = evt.getValidHandle<KalRepPtrCollection>(std::get<1>(vetoDefs_[iveto]));

      auto& cc = *std::get<2>(vetoDefs_[iveto]);
      int numVetoCandidates = 0;
      for(const auto& ptr: *uemh) {
        if(cc.accepted(ptr, evt, kdiag_, wh_)) {
          ++numVetoCandidates;
        }
      }

      hNumVetoCandidates_->Fill(iveto, numVetoCandidates, wh_.weight());
      if(numVetoCandidates>0) {
        return EventCutNumber::TrackVeto;
      }

    }

    return EventCutNumber::accepted;

  } // processTrack()

  //================================================================

} // namespace mu2e
