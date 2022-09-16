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
  struct ClusterState{
    ClusterState(Chisq const& chi2, WHSCOL const& hitstates) : chi2_(chi2), hitstates_(hitstates) {}
    ClusterState() {}
    Chisq chi2_; // total chisquared for this cluster
    WHSCOL hitstates_; // states for all hits in the cluster
    // merge with another score
    void merge(ClusterState const& other);
  };
  using ClusterStateCOL = std::vector<ClusterState>;

  class CombinatoricStrawHitUpdater {
    public:
      using CSHUConfig = std::tuple<unsigned,float,float,float,float,std::string,std::string,bool,int>;
      static std::string const& configDescription(); // description of the variables
      // struct to sort hit states by chisquared value
      struct ClusterStateComp {
        bool operator()(ClusterState const& a, ClusterState const& b)  const {
          return a.chi2_.chisqPerNDOF() < b.chi2_.chisqPerNDOF();
        }
      };
      CombinatoricStrawHitUpdater(CSHUConfig const& cshuconfig);
      ClusterState selectBest(ClusterStateCOL& cscores) const; // find the best cluster configuration given the score for each
      auto inactivePenalty() const { return inactivep_;}
      auto nullPenalty() const { return nullp_;}
      auto const& allowed() const { return allowed_; }
      auto minDeltaChi2() const { return mindchi2_; }
      auto nullDOCA() const { return nulldoca_; }
      // the work is done here
      template <class KTRAJ> void updateCluster(KKStrawHitCluster<KTRAJ>& cluster,KinKal::MetaIterConfig const& miconfig) const;
    private:
      unsigned csize_ =0; // minimum cluster size to update
      double inactivep_ =0; // chisquared penalty for inactive hits
      double nullp_ =0; // chisquared penalty for null hits
      double mindchi2_ =0; // minimum chisquared separation to consider 'significant'
      double nulldoca_ =0; // DOCA used to set null hit variance
      WHSMask freeze_; // states to freeze
      bool unfreeze_ =false; // ignore freeze state on input
      int diag_ =0; // diag print level
      WHSCOL allowed_; // allowed states
  };
  std::ostream& operator <<(std::ostream& os, ClusterState const& cscore );
}
#endif
