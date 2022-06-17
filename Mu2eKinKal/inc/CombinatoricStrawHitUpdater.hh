#ifndef Mu2eKinKal_CombinatoricStrawHitUpdater_hh
#define Mu2eKinKal_CombinatoricStrawHitUpdater_hh
//
//  StrawHitGroup updating using an exhaustive combinatoric algorithm, following the BTrk PanelAmbigResolver algorithm
//
#include "KinKal/Detector/WireHitStructs.hh"
#include "KinKal/General/Chisq.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/WHSIterator.hh"
#include <tuple>
#include <vector>
#include <memory>
#include <iostream>

namespace mu2e {
  using KinKal::WireHitState;
  using KinKal::Chisq;
  using KinKal::WireHitState;
  using KinKal::Parameters;
  using KinKal::Weights;
  using WHSCOL = std::vector<WireHitState>;
  using CHI2WHS = std::tuple<Chisq,WHSCOL>;
  using CHI2WHSCOL = std::vector<CHI2WHS>;

  class CombinatoricStrawHitUpdater {
    public:
      struct CHI2Comp { // comparator to sort hit states by chisquared value
        bool operator()(CHI2WHS const& a, CHI2WHS const& b)  const {
          return std::get<0>(a).chisqPerNDOF() < std::get<0>(b).chisqPerNDOF();
        }
      };
      CombinatoricStrawHitUpdater(double inactivep, double nullp, double mindchi2) :
        inactivep_(inactivep), nullp_(nullp), mindchi2_(mindchi2),
        allowed_{WireHitState::inactive, WireHitState::left, WireHitState::null, WireHitState::right} {}
      WHSCOL selectBest(CHI2WHSCOL& chi2s) const; // find the best configuration given the total chisq for each
      double penalty(WireHitState const& whs) const; // compute the penalty for each hit in a given state
      double inactivePenalty() const { return inactivep_;}
      double nullPenalty() const { return nullp_;}
      auto const& allowed() const { return allowed_; }
      double minDeltaChi2() const { return mindchi2_; }
// the work is done here
      template <class KTRAJ> void updateHits(std::vector<std::shared_ptr<KKStrawHit<KTRAJ>>>& hits) const;
    private:
      double inactivep_; // chisquared penalty for inactive hits
      double nullp_; // chisquared penalty for null hits
      double mindchi2_; // minimum chisquared separation to consider 'significant'
      WHSCOL allowed_; // allowed states
  };

  template<class KTRAJ> void CombinatoricStrawHitUpdater::updateHits(std::vector<std::shared_ptr<KKStrawHit<KTRAJ>>>& hits) const {
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
      uweights = Weights(uparams);
    }

    // iterate over all possible states of each hit, and compute the total chisquared for that
    WHSIterator whsiter(hits.size(),allowed());
    size_t ncombo = whsiter.nCombo();
    CHI2WHSCOL chi2whscol(0);
    chi2whscol.reserve(ncombo);
    do {
      // loop over hits in order and compute their residuals
      double chisq(0.0);
      unsigned ndof(0);
      for(size_t ihit=0;ihit < hits.size(); ++ihit) {
        auto const& shptr = hits[ihit];
        auto const& whstate = whsiter.current()[ihit];
        if(whstate != WireHitState::inactive) {
          // compute the chisquared contribution for this hit against the current parameters
          RESIDCOL resids;
          // compute residuals using this state (still WRT the reference parameters)
          DVEC dpvec = uparams.parameters() - shptr->referenceParameters().parameters();
          shptr->setResiduals(whstate,resids);
          for(auto resid : resids) {
            if(resid.active()){
              // update residuals to refer to unbiased parameters
              double uresidval = resid.value() - ROOT::Math::Dot(dpvec,resid.dRdP());
              double pvar = ROOT::Math::Similarity(resid.dRdP(),uparams.covariance());
              if(pvar<0) throw cet::exception("RECO")<<"mu2e::KKStrawHitGroup: negative variance " << std::endl;
              Residual uresid(uresidval,resid.variance(),pvar,resid.active(),resid.dRdP());
              chisq += uresid.chisq();
              ++ndof;
            }
          }
          // add null penalty
          if(whstate == WireHitState::null) chisq += nullp_;
          // compute the weight from this hits state, and update the parameters to use for subsequent hits TODO
        } else {
        // add penalty
          chisq += inactivep_;
          ++ndof; // count this as a DOF
        }
      }
      std::cout << "chisquared, ndof " << chisq << ndof << " States ";
      for(auto whstate : whsiter.current()) std::cout << "  " << whstate.state_;
      std::cout << std::endl;
      // record this chisq with the state of all the hits
      chi2whscol.emplace_back(Chisq(chisq,ndof),whsiter.current());
    } while(whsiter.increment());
    // test
    if(chi2whscol.size() != ncombo){
      throw cet::exception("RECO")<<"mu2e::KKStrawHitGroup: incomplete chisquared combinatorics" << std::endl;
    }
    // Let the updater choose the best state
    WHSCOL best = selectBest(chi2whscol);
    // assign the individual hit states according to this
    for(size_t ihit=0;ihit < hits.size(); ++ihit) {
      auto const& shptr = hits[ihit];
      auto const& whstate = best[ihit];
      shptr->setState(whstate);
    }
  }
}
#endif
