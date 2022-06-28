#ifndef Mu2eKinKal_CombinatoricStrawHitUpdater_hh
#define Mu2eKinKal_CombinatoricStrawHitUpdater_hh
//
//  StrawHitCluster updating using an exhaustive combinatoric algorithm, following the BTrk PanelAmbigResolver algorithm
//
#include "KinKal/General/Chisq.hh"
#include "KinKal/Fit/MetaIterConfig.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/WHSIterator.hh"
#include "Offline/Mu2eKinKal/inc/WireHitState.hh"
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
      CombinatoricStrawHitUpdater(double inactivep, double nullp, double mindchi2,int diag=0) :
        inactivep_(inactivep), nullp_(nullp), mindchi2_(mindchi2),diag_(diag),
        allowed_{WireHitState::inactive, WireHitState::left, WireHitState::null, WireHitState::right} {}
      ClusterScore selectBest(ClusterScoreCOL& cscores) const; // find the best cluster configuration given the score for each
      double penalty(WireHitState const& whs) const; // compute the penalty for each hit in a given state
      double inactivePenalty() const { return inactivep_;}
      double nullPenalty() const { return nullp_;}
      auto const& allowed() const { return allowed_; }
      double minDeltaChi2() const { return mindchi2_; }
      // the work is done here
      template <class KTRAJ> void updateHits(std::vector<std::shared_ptr<KKStrawHit<KTRAJ>>>& hits,KinKal::MetaIterConfig const& miconfig) const;
    private:
      double inactivep_; // chisquared penalty for inactive hits
      double nullp_; // chisquared penalty for null hits
      double mindchi2_; // minimum chisquared separation to consider 'significant'
      int diag_; // diag print level
      WHSCOL allowed_; // allowed states
      double wireHitRank(WHSCOL const& whscol) const; // rank wire hit states by 'conservativeness'
  };

  template<class KTRAJ> void CombinatoricStrawHitUpdater::updateHits(
      std::vector<std::shared_ptr<KKStrawHit<KTRAJ>>>& allhits,
      KinKal::MetaIterConfig const& miconfig) const {
    using SHCOL = std::vector<std::shared_ptr<KKStrawHit<KTRAJ>>>;
    // leave alone hits that were forced inactive
    SHCOL hits;
    hits.reserve(allhits.size());
    for(auto& shptr : allhits)
      if(shptr->hitState().usable())
        hits.push_back(shptr);
    //
    if(hits.size() == 0)return;
    // sort the hit ptrs by time
    std::sort(hits.begin(),hits.end(),StrawHitTimeSort<KTRAJ>());
    Parameters uparams;
    Weights uweights;
    // Find the first active hit
    auto early = hits.end();
    for(auto ish=hits.begin(); ish != hits.end(); ++ish){
      if((*ish)->active()){
        early = ish;
        break;
      }
    }
    if(early != hits.end()){
      uweights = Weights((*early)->referenceParameters());
      for (auto const& sh : hits) {
        if(sh->active()) uweights -= sh->weight();
      }
      uparams = Parameters(uweights);
    } else {
      uparams = hits.front()->referenceParameters();
    }
    uweights = Weights(uparams);
   // check determinant
    double determinant;
    if(!uparams.covariance().Det(determinant) || determinant < std::numeric_limits<float>::min()){
      if(diag_ > 0)std::cout << "Negative unbiased covar determinant = " << determinant << std::endl;
// for now, nullify these hits and move on
// for now, do nothing
//      for(size_t ihit=0;ihit < hits.size(); ++ihit) {
//        static WireHitState nullstate(WireHitState::null);
//        auto const& shptr = hits[ihit];
//        shptr->setState(nullstate);
//      }
//      TODO: try a more careful weight subtraction including material effects.
      return;
    }

    // iterate over all possible states of each hit, and compute the total chisquared for that
    WHSIterator whsiter(hits.size(),allowed());
    size_t ncombo = whsiter.nCombo();
    ClusterScoreCOL cscores(0);
    cscores.reserve(ncombo);
    do {
      auto cparams = uparams;
      auto cweights = uweights;
      // loop over hits in order and compute their residuals
      double chisq(0.0);
      unsigned ndof(0);
      for(size_t ihit=0;ihit < hits.size(); ++ihit) {
        auto const& shptr = hits[ihit];
        auto const& whstate = whsiter.current()[ihit];
        if(whstate.active()) {
          // compute the chisquared contribution for this hit against the current parameters
          RESIDCOL resids;
          // compute residuals using this state (still WRT the reference parameters)
          DVEC dpvec = cparams.parameters() - shptr->referenceParameters().parameters();
          shptr->setResiduals(whstate,resids);
          for(auto resid : resids) {
            if(resid.active()){
              // update residuals to refer to unbiased parameters
              double uresidval = resid.value() - ROOT::Math::Dot(dpvec,resid.dRdP());
              double pvar = ROOT::Math::Similarity(resid.dRdP(),cparams.covariance());
//              if(pvar<0) throw cet::exception("RECO")<<"mu2e::KKStrawHitCluster: negative variance " << pvar << std::endl;
              if(pvar<0) std::cout <<"mu2e::KKStrawHitCluster: negative variance " << pvar
                << " determinant = " << determinant << std::endl;
              Residual uresid(uresidval,resid.variance(),pvar,resid.active(),resid.dRdP());
              chisq += uresid.chisq();
              ++ndof;
            }
          }
          // add null penalty
          if(whstate == WireHitState::null) chisq += nullp_;
          // compute the weight from this hits state, and update the parameters to use for subsequent hits
          // this isn't necessary for the last hit since there are no subsequent hits
          if(ihit+1 < hits.size()){
            for(auto resid : resids) {
              if(resid.active())cweights +=  resid.weight(cparams.parameters(),miconfig.varianceScale());
            }
            cparams = Parameters(cweights);
          }
        } else {
        // add penalty
          chisq += inactivep_;
          ++ndof; // count this as a DOF
        }
      }
      if(diag_ > 1){
        std::cout << "chi2 " << chisq << " ndof " << ndof << " States ";
        for(auto whstate : whsiter.current()) std::cout << "  " << whstate.state_;
        std::cout << std::endl;
      }
      // record this chisq with the state of all the hits
      cscores.emplace_back(Chisq(chisq,ndof),whsiter.current());
    } while(whsiter.increment());
    // test
    if(cscores.size() != ncombo){
      throw cet::exception("RECO")<<"mu2e::KKStrawHitCluster: incomplete chisquared combinatorics" << std::endl;
    }
    // choose the best cluster score
    auto best = selectBest(cscores);
    if(diag_ > 0){
      std::cout << "Best Combo chi2 " << best.chi2_ << " states ";
      for(auto whs : best.hitstates_)std::cout << "  " << whs.state_;
      std::cout << std::endl;
    }
    // assign the individual hit states according to this
    for(size_t ihit=0;ihit < hits.size(); ++ihit) {
      hits[ihit]->setState(best.hitstates_[ihit]);
    }
  }
}
#endif
