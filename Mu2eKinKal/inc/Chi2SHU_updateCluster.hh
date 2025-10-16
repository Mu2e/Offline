//
//  Implementation of the hit updater function for Chi2SHU
//
#ifndef Mu2eKinKal_Chi2SHU_updateCluster_hh
#define Mu2eKinKal_Chi2SHU_updateCluster_hh
#include "Offline/Mu2eKinKal/inc/Chi2SHU.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHitCluster.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
namespace mu2e {
  template<class KTRAJ> void Chi2SHU::updateCluster(
      KKStrawHitCluster<KTRAJ>& shcluster,
      KinKal::MetaIterConfig const& miconfig) const {
    using KKSTRAWHITCLUSTER = KKStrawHitCluster<KTRAJ>;
    using KKSTRAWHIT = KKStrawHit<KTRAJ>;
    using SHCOL = std::vector<std::shared_ptr<KKSTRAWHIT>>;
    // only act on clusters with a minimum # of hits
    // dont act on hits that are unusable or frozen
    SHCOL hits;
    hits.reserve(shcluster.strawHits().size());
    for(auto& shptr : shcluster.strawHits()){
      if(shptr->hitState().updateable(StrawHitUpdaters::Chi2) || ( shptr->hitState().usable() && shptr->hitState().isIn(unfreeze_))) hits.push_back(shptr);
    }
    // make sure this cluster meets the requirements for updating
    if(hits.size() < csize_ )return;
    // get the reference weight as starting point for the unbiased weight + parameters
    Weights uweights = Weights(hits.front()->referenceParameters());
    // subtract the weight of active hits from this reference; this removes their bias from the reference
    for (auto const& sh : hits) {
      if(sh->active()) uweights -= sh->weight();
    }
    // invert to get unbiased parameters
    Parameters uparams = Parameters(uweights);
    // in marginal fits removing even a few hits can leave the fit underconstrained, resulting in a zero determinant.  For now, do nothing
    // with these clusters
    double determinant;
    if(!uparams.covariance().Det(determinant) || determinant < std::numeric_limits<float>::min()){
      if(diag_ > 2)std::cout << "Negative unbiased covar determinant = " << determinant << std::endl;
      return;
    }
    // iterate over allowed states of each hit, and incrementally compute the total chisquared for all the hits in the cluster WRT the unbiased parameters
    WHSIterator whsiter(hits.size(),allowed_);
    // set the null hit docas: it can't be done in the
    size_t nstates = whsiter.nStates();
    ClusterStateCOL cstates(0);
    cstates.reserve(nstates);
    do {
      auto cparams = uparams;
      auto cweights = uweights;
      // loop over hits and compute their residuals
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
          // only use distance residual
          auto const& resid = resids[Mu2eKinKal::dresid];
          if(resid.active()) {
            // update residuals to refer to unbiased parameters
            double uresidval = resid.value() - ROOT::Math::Dot(dpvec,resid.dRdP());
            double pvar = ROOT::Math::Similarity(resid.dRdP(),cparams.covariance());
            //              if(pvar<0) throw cet::exception("RECO")<<"mu2e::KKStrawHitCluster: negative variance " << pvar << std::endl;
            if(pvar<0){
              // another symptom of under-constrained clusters is negative variances.  I need a better strategy for these TODO
              if(diag_ > 2) std::cout <<"mu2e::KKStrawHitCluster: negative variance " << pvar
                << " determinant = " << determinant << std::endl;
              pvar = resid.parameterVariance();
            }
            // Use the unbiased residual to compute the chisq
            Residual uresid(uresidval,resid.variance(),pvar,resid.dRdP(),resid.active());
            chisq += uresid.chisq();
            ++ndof;
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
      // record this chisq with the state of all the hits
      cstates.emplace_back(Chisq(chisq,ndof),whsiter.current());
    } while(whsiter.increment());
    // test
    if(cstates.size() != nstates){
      throw cet::exception("RECO")<<"mu2e::KKStrawHitCluster: incomplete chisquared combinatorics" << std::endl;
    }
    // choose the best cluster state
    auto best = selectBest(cstates);
    // assign the individual hit states according to this, and update their fit info
    for(size_t ihit=0;ihit < hits.size(); ++ihit) {
      auto& whs = best.hitstates_[ihit];
      hits[ihit]->setState(whs);
    }
  }
}
#endif
