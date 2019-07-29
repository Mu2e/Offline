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

#ifndef CutAndCountAnalysis_inc_CutAndCountAnalysis_hh
#define CutAndCountAnalysis_inc_CutAndCountAnalysis_hh

#include <string>
#include <vector>

#include "art_root_io/TFileDirectory.h"

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

#include "TH1.h"
#include "TH2.h"

#include "CutAndCountAnalysis/inc/TrackLevelCuts.hh"

namespace mu2e {

  namespace CutAndCount {
    enum class EventCutNumber;

    // name, input collection, cuts
    typedef std::tuple<std::string,art::InputTag,std::unique_ptr<TrackLevelCuts> > TrackVetoDef;
  }

  class CutAndCountAnalysis {
  public:
    explicit CutAndCountAnalysis(const fhicl::ParameterSet& pset, art::TFileDirectory tfdir);

    bool accepted(const art::Event& event);

  private:
    //----------------------------------------------------------------
    art::InputTag signalCandidateInput_;
    TrackLevelCuts signalTrackCuts_;

    std::vector<CutAndCount::TrackVetoDef> vetoDefs_;

    EventWeightHelper wh_;
    KalDiag kdiag_;

    TH1 *h_cuts_p_;
    TH1 *h_cuts_r_;
    TH1 *w_cuts_p_;
    TH1 *w_cuts_r_;
    TH1 *hNumSignalCandidates_;
    TH2 *hNumVetoCandidates_;

    CutAndCount::EventCutNumber processEvent(const art::Event& evt);
  };

} // namespace mu2e

#endif/*CutAndCountAnalysis_inc_CutAndCountAnalysis_hh*/
