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
    // merge with another score
    void merge(ClusterScore const& other);
  };
  using ClusterScoreCOL = std::vector<ClusterScore>;

  class CombinatoricStrawHitUpdater {
    public:
      using CSHUConfig = std::tuple<unsigned,float,float,float,float,bool,bool,bool,int>;
      // struct to sort hit states by chisquared value
      struct ClusterScoreComp {
        bool operator()(ClusterScore const& a, ClusterScore const& b)  const {
          return a.chi2_.chisqPerNDOF() < b.chi2_.chisqPerNDOF();
        }
      };
      CombinatoricStrawHitUpdater(CSHUConfig const& cshuconfig) {
        csize_ = std::get<0>(cshuconfig);
        inactivep_ = std::get<1>(cshuconfig);
        nullp_ = std::get<2>(cshuconfig);
        mindchi2_ = std::get<3>(cshuconfig);
        nulldoca_ = std::get<4>(cshuconfig);
        allownull_ = std::get<5>(cshuconfig);
        nulltime_ = std::get<6>(cshuconfig);
        freeze_ = std::get<7>(cshuconfig);
        diag_ = std::get<8>(cshuconfig);
        static double invthree(1.0/3.0);
        dvar_ = invthree*nulldoca_*nulldoca_;
        if(allownull_)
          allowed_ = WHSCOL{WireHitState::inactive, WireHitState::left, WireHitState::null, WireHitState::right};
        else
          allowed_ = WHSCOL{WireHitState::inactive, WireHitState::left, WireHitState::right};
        // set the algorithm; this propagates to the StrawHits
        for(auto& whs : allowed_)whs.algo_ = StrawHitUpdaters::Combinatoric;
        toff_ = tvar_ = 0.0;
        std::cout << "CombinatoricStrawHitUpdater " << inactivep_ << " " << nullp_ << " " << mindchi2_ << " " << nulldoca_ << " " << allownull_ << " " << nulltime_ << std::endl;
      }
      ClusterScore selectBest(ClusterScoreCOL& cscores) const; // find the best cluster configuration given the score for each
      double penalty(WireHitState const& whs) const; // compute the penalty for each hit in a given state
      auto inactivePenalty() const { return inactivep_;}
      auto nullPenalty() const { return nullp_;}
      auto const& allowed() const { return allowed_; }
      auto minDeltaChi2() const { return mindchi2_; }
      auto nullDOCA() const { return nulldoca_; }
      auto allowNull() const { return allownull_; }
      auto nullTime() const { return nulltime_; }
      auto distVariance() const { return dvar_; }
      // the work is done here
      template <class KTRAJ> void updateCluster(KKStrawHitCluster<KTRAJ>& cluster,KinKal::MetaIterConfig const& miconfig) const;
    private:
      unsigned csize_; // minimum cluster size to update
      double inactivep_; // chisquared penalty for inactive hits
      double nullp_; // chisquared penalty for null hits
      double mindchi2_; // minimum chisquared separation to consider 'significant'
      double nulldoca_; // DOCA used to set null hit variance
      bool allownull_; // allow null ambiguity assignment
      bool nulltime_; // use time residual in null hit chisquared
      bool freeze_; // freeze disambiguated clusters
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
