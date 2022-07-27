#ifndef Mu2eKinKal_CombinatoricStrawHitUpdater_hh
#define Mu2eKinKal_CombinatoricStrawHitUpdater_hh
//
//  StrawHitCluster updating using an exhaustive combinatoric algorithm, following the BTrk PanelAmbigResolver algorithm
//
#include "KinKal/General/Chisq.hh"
#include "KinKal/Fit/MetaIterConfig.hh"
#include "KinKal/General/Parameters.hh"
#include "KinKal/General/Weights.hh"
#include "Offline/Mu2eKinKal/inc/WHSIterator.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdater.hh"
#include <tuple>
#include <vector>
#include <memory>
#include <limits>
#include <iostream>

namespace mu2e {
  using KinKal::Chisq;
  using KinKal::Parameters;
  using KinKal::Weights;
  using WHSCOL = std::vector<WireHitState>;
  template<class KTRAJ> class KKStrawHitCluster;
  struct ClusterScore{
    ClusterScore(Chisq const& chi2, WHSCOL const& hitstates) : chi2_(chi2), hitstates_(hitstates) {}
    ClusterScore() {}
    Chisq chi2_; // total chisquared for this cluster
    WHSCOL hitstates_; // states for all hits in the cluster
  };
  using ClusterScoreCOL = std::vector<ClusterScore>;

  class CombinatoricStrawHitUpdater : public StrawHitUpdater {
    public:
      using CSHUConfig = std::tuple<float,float,float,float,bool,bool,int>;
      // struct to sort hit states by chisquared value
      struct ClusterScoreComp {
        bool operator()(ClusterScore const& a, ClusterScore const& b)  const {
          return a.chi2_.chisqPerNDOF() < b.chi2_.chisqPerNDOF();
        }
      };
      CombinatoricStrawHitUpdater(CSHUConfig const& cshuconfig) {
        inactivep_ = std::get<0>(cshuconfig);
        nullp_ = std::get<1>(cshuconfig);
        mindchi2_ = std::get<2>(cshuconfig);
        mindoca_ = std::get<3>(cshuconfig);
        allownull_ = std::get<4>(cshuconfig);
        nulltime_ = std::get<5>(cshuconfig);
        diag_ = std::get<6>(cshuconfig);
        static double invthree(1.0/3.0);
        dvar_ = invthree*mindoca_*mindoca_;
        if(allownull_)
          allowed_ = WHSCOL{WireHitState::inactive, WireHitState::left, WireHitState::null, WireHitState::right};
        else
          allowed_ = WHSCOL{WireHitState::inactive, WireHitState::left, WireHitState::right};
        toff_ = tvar_ = 0.0;
        std::cout << "CombinatoricStrawHitUpdater " << inactivep_ << " " << nullp_ << " " << mindchi2_ << " " << mindoca_ << " " << allownull_ << " " << nulltime_ << std::endl;
      }
      // base class interface: note the WireHitState function isn't used
      WireHitState wireHitState(ClosestApproachData const& tpdata,Straw const& straw) const override { return WireHitState(WireHitState::inactive); }
      NullHitInfo nullHitInfo(StrawResponse const& sresponse,Straw const& straw) const override;
      StrawHitUpdaters::algorithm algorithm() const override { return StrawHitUpdaters::Combinatoric; }
      bool useUnbiasedClosestApproach() const override { return false; }
      bool useStrawHitCluster() const override { return true; }

      ClusterScore selectBest(ClusterScoreCOL& cscores) const; // find the best cluster configuration given the score for each
      double penalty(WireHitState const& whs) const; // compute the penalty for each hit in a given state
      auto inactivePenalty() const { return inactivep_;}
      auto nullPenalty() const { return nullp_;}
      auto const& allowed() const { return allowed_; }
      auto minDeltaChi2() const { return mindchi2_; }
      auto minDOCA() const { return mindoca_; }
      auto allowNull() const { return allownull_; }
      auto nullTime() const { return nulltime_; }
      auto distVariance() const { return dvar_; }
      // the work is done here
      template <class KTRAJ> void updateCluster(KKStrawHitCluster<KTRAJ>& cluster,KinKal::MetaIterConfig const& miconfig) const;
    private:
      double inactivep_; // chisquared penalty for inactive hits
      double nullp_; // chisquared penalty for null hits
      double mindchi2_; // minimum chisquared separation to consider 'significant'
      double mindoca_; // minimum DOCA for LR ambiguity assigned hits
      bool allownull_; // allow null ambiguity assignment
      bool nulltime_; // use time residual in null hit chisquared
      int diag_; // diag print level
      WHSCOL allowed_; // allowed states
      double wireHitRank(WHSCOL const& whscol) const; // rank wire hit states by 'conservativeness'
      double dvar_; // null hit distance variance
      mutable double toff_; // null hit time offset cache
      mutable double tvar_; // null hit time variance cache
  };
  std::ostream& operator <<(std::ostream& os, ClusterScore const& cscore );
}
#endif
