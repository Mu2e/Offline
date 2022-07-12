//
//  Implementation of the hit updater function for CombinatoricStrawHitUpdater
//
#ifndef Mu2eKinKal_CSHU_updateCluster_hh
#define Mu2eKinKal_CSHU_updateCluster_hh
#include "Offline/Mu2eKinKal/inc/CombinatoricStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
namespace mu2e {
  template<class KTRAJ> void CombinatoricStrawHitUpdater::updateCluster(
      KKStrawHitCluster<KTRAJ>& shcluster,
      KinKal::MetaIterConfig const& miconfig) const {
    using KKSTRAWHITCLUSTER = KKStrawHitCluster<KTRAJ>;
    using KKSTRAWHIT = KKStrawHit<KTRAJ>;
    using SHCOL = std::vector<std::shared_ptr<KKSTRAWHIT>>;
    // leave alone hits that were forced inactive
    SHCOL hits;
    hits.reserve(shcluster.strawHits().size());
    for(auto& shptr : shcluster.strawHits())
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
      // subtract the material; this requires obtaining the covariance from the
      // xing, which needs some restructuring TODO
    } else {
      uparams = hits.front()->referenceParameters();
    }
    uweights = Weights(uparams);
   // check determinant
    double determinant;
    if(!uparams.covariance().Det(determinant) || determinant < std::numeric_limits<float>::min()){
      if(diag_ > 0)std::cout << "Negative unbiased covar determinant = " << determinant << std::endl;
      // without these hits the fit is underconstrained so the algorithm fails
      // for now, do nothing
      //      for(size_t ihit=0;ihit < hits.size(); ++ihit) {
      //        static WireHitState nullstate(WireHitState::null);
      //        auto const& shptr = hits[ihit];
      //        shptr->setState(nullstate);
      //      }
      //      TODO: try a more careful weight subtraction including material effects.
      //      Or, revert to a hit-by-hit update of this cluster
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
          shptr->setResiduals(miconfig,whstate,resids);
          for(size_t iresid= 0; iresid < resids.size(); ++iresid) {
            auto const& resid = resids[iresid];
            if(resid.active() && (whstate.useDrift() || (iresid==KKSTRAWHIT::dresid || nulltime_))){
             // update residuals to refer to unbiased parameters
              double uresidval = resid.value() - ROOT::Math::Dot(dpvec,resid.dRdP());
              double pvar = fabs(ROOT::Math::Similarity(resid.dRdP(),cparams.covariance()));
//              if(pvar<0) throw cet::exception("RECO")<<"mu2e::KKStrawHitCluster: negative variance " << pvar << std::endl;
              if(diag_ > 0 && pvar<0) std::cout <<"mu2e::KKStrawHitCluster: negative variance " << pvar
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
    // assign the individual hit states according to this, and update their fit info
    for(size_t ihit=0;ihit < hits.size(); ++ihit) {
      hits[ihit]->setState(miconfig,best.hitstates_[ihit]);
    }
  }
}
#endif
