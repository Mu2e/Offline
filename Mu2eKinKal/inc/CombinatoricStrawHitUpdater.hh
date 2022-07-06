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
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
#include "Offline/Mu2eKinKal/inc/StrawHitUpdaters.hh"
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
  template<class KTRAJ> class KKStrawHit;
  struct ClusterScore{
    ClusterScore(Chisq const& chi2, WHSCOL const& hitstates) : chi2_(chi2), hitstates_(hitstates) {}
    ClusterScore() {}
    Chisq chi2_; // total chisquared for this cluster
    WHSCOL hitstates_; // states for all hits in the cluster
  };
  using ClusterScoreCOL = std::vector<ClusterScore>;

  class CombinatoricStrawHitUpdater {
    public:
      // struct to sort hit states by chisquared value
      struct ClusterScoreComp {
        bool operator()(ClusterScore const& a, ClusterScore const& b)  const {
          return a.chi2_.chisqPerNDOF() < b.chi2_.chisqPerNDOF();
        }
      };
      CombinatoricStrawHitUpdater(double inactivep, double nullp, double mindchi2,double nulldoca,int diag=0) :
        inactivep_(inactivep), nullp_(nullp), mindchi2_(mindchi2), nulldoca_(nulldoca), diag_(diag),
        allowed_{WireHitState::inactive, WireHitState::left, WireHitState::null, WireHitState::right} {}
      ClusterScore selectBest(ClusterScoreCOL& cscores) const; // find the best cluster configuration given the score for each
      double penalty(WireHitState const& whs) const; // compute the penalty for each hit in a given state
      auto inactivePenalty() const { return inactivep_;}
      auto nullPenalty() const { return nullp_;}
      auto const& allowed() const { return allowed_; }
      auto minDeltaChi2() const { return mindchi2_; }
      auto meanNullDOCA() const { return nulldoca_; }
      StrawHitUpdaters::algorithm algorithm() const { return StrawHitUpdaters::Combinatoric; }
      // the work is done here
      template <class KTRAJ> void updateHits(std::vector<std::shared_ptr<KKStrawHit<KTRAJ>>>& hits,KinKal::MetaIterConfig const& miconfig) const;
    private:
      double inactivep_; // chisquared penalty for inactive hits
      double nullp_; // chisquared penalty for null hits
      double mindchi2_; // minimum chisquared separation to consider 'significant'
      double nulldoca_; // effective mean DOCA for null ambiguity assigned hits
      int diag_; // diag print level
      WHSCOL allowed_; // allowed states
      double wireHitRank(WHSCOL const& whscol) const; // rank wire hit states by 'conservativeness'
  };

}
#endif
