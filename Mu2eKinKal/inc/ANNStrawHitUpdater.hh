//
// ANN-based  ambiguity updater.  This assigns LR ambiguity based on the predicted accuracy of the DOCA sign
//
#ifndef Mu2eKinKal_ANNStrawHitUpdater_hh
#define Mu2eKinKal_ANNStrawHitUpdater_hh
#include "Offline/Mu2eUtilities/inc/MVATools.hh"
#include "KinKal/Trajectory/ClosestApproachData.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/WHSMask.hh"
#include "Offline/Mu2eKinKal/inc/DriftInfo.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
#include <tuple>
#include <string>
#include <iostream>

namespace mu2e {
  class ComboHit;
  // Update based just on ANN to the wire
  class ANNStrawHitUpdater {
    public:
      using ANNSHUConfig = std::tuple<std::string,float,float,std::string>;
      static std::string const& configDescription(); // description of the variables
      ANNStrawHitUpdater() : mva_(0), mvacut_(0.0), nulldoca_(2.5) {}
      ANNStrawHitUpdater(ANNStrawHitUpdater const& other) : mva_(0), mvacut_(other.mvacut_), nulldoca_(other.nulldoca_), freeze_(other.freeze_),
      mintdrift_(other.mintdrift_), maxtdrift_(other.maxtdrift_), maxdoca_(other.maxdoca_), maxresidpull_(other.maxresidpull_) {
        if(other.mva_) mva_ = new MVATools(*other.mva_);
      }
      ~ANNStrawHitUpdater() { delete mva_; }
      ANNStrawHitUpdater(ANNSHUConfig const& annshuconfig);
      WireHitState wireHitState(WireHitState const& input, KinKal::ClosestApproachData const& tpdata, DriftInfo const& dinfo, ComboHit const& chit) const;
    private:
      MVATools* mva_; // neural net calculator
      double mvacut_; // cut value to decide if drift information is usable
      double nulldoca_; // null hit doca
      WHSMask freeze_; // states to freeze
      double mintdrift_,maxtdrift_, maxdoca_; // outlier cuts
      double maxresidpull_;
  };
}
#endif
