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

#include "art/Framework/Services/Optional/TFileDirectory.h"

#include "fhiclcpp/ParameterSet.h"
#include "art/Framework/Principal/Event.h"
#include "art/Utilities/InputTag.h"

#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "BTrk/KalmanTrack/KalRep.hh"
#include "TrkDiag/inc/KalDiag.hh"

#include "Mu2eUtilities/inc/EventWeightHelper.hh"

#include "TrackCaloMatching/inc/TrackClusterMatch.hh"
//#include "RecoDataProducts/inc/TrkCaloMatchCollection.hh"
//#include "RecoDataProducts/inc/TrkCaloMatch.hh"
#include "RecoDataProducts/inc/TrkCaloIntersect.hh"

#include "TH1.h"
#include "TH2.h"

#include "CutAndCountAnalysis/inc/TrkCutHist.hh"

namespace mu2e {

  namespace CutAndCount {
    enum class TrkCut;
  }

  class CutAndCountAnalysis {
  public:
    explicit CutAndCountAnalysis(const fhicl::ParameterSet& pset, art::TFileDirectory tfdir);

    bool accepted(const art::Event& event);

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

      // track-calo match
      double caloMatchChi2;
      double caloemin;
      double caloemax;

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

        , caloMatchChi2(pset.get<double>("caloMatchChi2"))
        , caloemin(pset.get<double>("caloemin"))
        , caloemax(pset.get<double>("caloemax"))

        , pmin(pset.get<double>("pmin"))
        , pmax(pset.get<double>("pmax"))
      {}
    };

    //----------------------------------------------------------------
    typedef std::vector<art::InputTag> InputTags;
    art::InputTag trackDemInput_;
    art::InputTag caloMatchDemInput_;
    PhysicsCuts cuts_;
    EventWeightHelper wh_;
    KalDiag kdiag_;

    TH1 *h_cuts_p_;
    TH1 *h_cuts_r_;
    TH1 *w_cuts_p_;
    TH1 *w_cuts_r_;
    TH1 *hNumAcceptedTracks_;
    CutAndCount::TrkCutHist hTrkCuts_;

    CutAndCount::TrkCut processTrack(const art::Ptr<KalRep>& trk, const art::Event& evt);
    const TrackClusterMatch* findCaloMatch(const art::Ptr<KalRep>& trk, const art::Event& evt);
  };

} // namespace mu2e

#endif/*CutAndCountAnalysis_inc_CutAndCountAnalysis_hh*/
