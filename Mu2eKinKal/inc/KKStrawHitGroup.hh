#ifndef Mu2eKinKal_KKStrawHitGroup_hh
#define Mu2eKinKal_KKStrawHitGroup_hh
//
//  class describing a group of nearby (in time and space) straw hits which need to be updated coherently.
//  Even though it is formally a 'Hit', its only interaction with the fit is in iteration updating.
//
#include "KinKal/Detector/Hit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/CombinatoricStrawHitUpdater.hh"
#include "Offline/Mu2eKinKal/inc/WHSIterator.hh"
#include "Offline/DataProducts/inc/StrawIdMask.hh"
#include "cetlib_except/exception.h"
#include <vector>
#include <memory>
#include <limits>
namespace mu2e {

  //
  // Functor to define which StrawHits belong in a single group
  //
  template <class KTRAJ> class KKStrawHitGrouper {
    public:
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      KKStrawHitGrouper( StrawIdMask groupLevel, double groupDt) :
        groupLevel_(groupLevel), groupDt_(groupDt) {}
      // decide if 2 hits should be grouped
      bool group(KKSTRAWHIT const& sh1, KKSTRAWHIT const& sh2) const {
        return groupLevel_.equal(sh1.strawId(),sh2.strawId()) &&
          fabs(sh1.time()-sh2.time()) < groupDt_;
      }
      auto groupLevel() const { return groupLevel_.level(); }
      double groupDt() const { return groupDt_; }
    private:
      StrawIdMask groupLevel_; // level to group straw hits
      double groupDt_; // max particle time difference between straw hits in a group
  };

  template <class KTRAJ> class KKStrawHitGroup : public KinKal::Hit<KTRAJ> {
    public:
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
      using KKSTRAWHITCOL = std::vector<KKSTRAWHITPTR>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using KKSTRAWHITGROUPER = KKStrawHitGrouper<KTRAJ>;
      // sort hits by time
      struct StrawHitSort {
        bool operator ()( const KKSTRAWHITPTR& hit1, const KKSTRAWHITPTR& hit2) {
          return hit1->time() < hit2->time(); }
      };
      KKStrawHitGroup() {}
      // create from a single hit
      KKStrawHitGroup(KKSTRAWHITPTR const& hitptr);
      // create from a collection of panel hits
      KKStrawHitGroup(KKSTRAWHITCOL const& hits,KKSTRAWHITGROUPER const& grouper);
      //Hit interface
      bool active() const override { return false; } // panel hits are never active
      KinKal::Chisq chisq(KinKal::Parameters const& params) const override { return KinKal::Chisq(); }
      unsigned nDOF() const override { return 0; }
      KinKal::Weights const& weight() const override { return (*hits_.begin())->weight(); }
      double time() const override;
      void updateReference(KTRAJPTR const& ktrajptr) override {} // nothing to do here, ref comes from individual hits
      KTRAJPTR const& refTrajPtr() const override { return (*hits_.begin())->refTrajPtr(); }
      // update the internals of the hit, specific to this meta-iteraion.  This will affect the next fit iteration
      void updateState(KinKal::MetaIterConfig const& config,bool first) override;
      void updateWeight(MetaIterConfig const& config) override {} // to be removed FIXME
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      ~KKStrawHitGroup(){}
      // KKStrawHitGroup specific interface
      auto const& strawHits() const { return hits_; }
      bool canAddHit(KKSTRAWHITPTR hit,KKSTRAWHITGROUPER const& grouper) const;
      void addHit(KKSTRAWHITPTR hit,KKSTRAWHITGROUPER const& grouper);
      // updater functions
      void updateCombo(const CombinatoricStrawHitUpdater* cshu);
    private:
      // references to the individual hits in this hit group
      KKSTRAWHITCOL hits_;
  };

  template<class KTRAJ>  KKStrawHitGroup<KTRAJ>::KKStrawHitGroup(KKSTRAWHITPTR const& hitptr) {
    hits_.push_back(hitptr);
  }

  template<class KTRAJ>  KKStrawHitGroup<KTRAJ>::KKStrawHitGroup(KKSTRAWHITCOL const& hits,KKSTRAWHITGROUPER const& grouper) : hits_(hits) {
    // verify grouping
    for (auto ihit=hits_.begin(); ihit != hits_.end(); ++ihit){
      for (auto jhit=ihit+1; jhit != hits_.end(); ++jhit){
        if(!grouper.group(**ihit,**jhit))
          throw cet::exception("RECO")<<"mu2e::KKStrawHitGroup: KKStrawHits don't belong in group"<< std::endl;
      }
    }
  }

  template<class KTRAJ> bool KKStrawHitGroup<KTRAJ>::canAddHit(KKSTRAWHITPTR hit, KKSTRAWHITGROUPER const& grouper) const {
    bool add(true);
    for(auto const& myhit : hits_)
      add &= grouper.group(*myhit,*hit);
    return add;
  }

  template<class KTRAJ> void KKStrawHitGroup<KTRAJ>::addHit(KKSTRAWHITPTR hit, KKSTRAWHITGROUPER const& grouper) {
    bool add = canAddHit(hit,grouper);
    if(!add)throw cet::exception("RECO")<<"mu2e::KKStrawHitGroup: KKStrawHit doesn't belong in group"<< std::endl;
    hits_.push_back(hit);
  }

  template<class KTRAJ> double KKStrawHitGroup<KTRAJ>::time() const {
    // return time just past the last hit's time.  This insures hit groups are updated after individual hits
    // that shouldn't matter, but...
    double maxtime(-std::numeric_limits<float>::max());
    unsigned nactive(0);
    for(auto const& hit : hits_){
      if(hit && hit->active())++nactive;
      maxtime = std::max(hit->time(),maxtime);
    }
    static double epsilon(1e-6);
    return maxtime + epsilon;
  }

  template<class KTRAJ> void KKStrawHitGroup<KTRAJ>::updateState(KinKal::MetaIterConfig const& miconfig,bool first) {
    if(first){
      // sort the hit ptrs
      std::sort(hits_.begin(),hits_.end(),StrawHitSort ());
      // look for an updater; if it's there, update the state
      // Extend this logic if new StrawHitGroup updaters are introduced
      auto cshu = miconfig.findUpdater<CombinatoricStrawHitUpdater>();
      if(cshu != 0){
        updateCombo(cshu);
      }
    }
  }

  // strawhit-specific functions assocated with Combinatoric updater.  These must go here as the Updater itself
  // can't be templated
  template<class KTRAJ> void KKStrawHitGroup<KTRAJ>::updateCombo(const CombinatoricStrawHitUpdater* cshu) {
    using KinKal::WireHitState;
    using KinKal::Parameters;
    using KinKal::Weights;
    Parameters uparams;
    Weights uweights;
    // Find the first active hit
    auto early = hits_.end();
    for(auto ish=hits_.begin(); ish != hits_.end(); ++ish){
      if((*ish)->active()){
        early = ish;
        break;
      }
    }
    if(early != hits_.end()){
      uweights = Weights((*early)->referenceParameters());
      for (auto const& sh : hits_) {
        if(sh->active()) uweights -= sh->weight();
      }
      uparams = Parameters(uweights);
    } else {
      uparams = hits_.front()->referenceParameters();
      uweights = Weights(uparams);
    }

    // iterate over all possible states of each hit, and compute the total chisquared for that
    WHSIterator whsiter(hits_.size(),cshu->allowed());
    size_t ncombo = whsiter.nCombo();
    CHI2WHSCOL chi2whscol(0);
    chi2whscol.reserve(ncombo);
    do {
      // loop over hits in order and compute their residuals
      double chisq(0.0);
      for(size_t ihit=0;ihit < hits_.size(); ++ihit) {
        auto const& shptr = hits_[ihit];
        auto const& whstate = whsiter.current()[ihit];
        if(whstate != WireHitState::inactive) {
          // compute the chisquared contribution for this hit against the current parameters
          RESIDCOL resids;
          // compute residuals using this state (still WRT the reference parameters)
          shptr->setResiduals(whstate,resids);
          for(auto resid : resids) {
            // update residuals to refer to unbiased parameters
            DVEC dpvec = uparams.parameters() - shptr->referenceParameters().parameters();
            double uresidval = resid.value() - ROOT::Math::Dot(dpvec,resid.dRdP());
            double pvar = ROOT::Math::Similarity(resid.dRdP(),uparams.covariance());
            if(pvar<0) throw cet::exception("RECO")<<"mu2e::KKStrawHitGroup: negative variance " << std::endl;
            Residual uresid(uresidval,resid.variance(),pvar,resid.active(),resid.dRdP());
            chisq += uresid.chisq();
          }
          // compute the weight from this hits state, and update the parameters for the next hit TODO
        }
        // add penalty terms
        chisq += cshu->penalty(whstate);
      }
      // record this chisq with the state of all the hits
      chi2whscol.emplace_back(chisq,whsiter.current());
    } while(whsiter.increment());
    // test
    if(chi2whscol.size() != ncombo){
     throw cet::exception("RECO")<<"mu2e::KKStrawHitGroup: incomplete chisquared combinatorics" << std::endl;
    }
    // Let the updater choose the best state
    WHSCOL best = cshu->selectBest(chi2whscol);
    // assign the individual hit states according to this
    for(size_t ihit=0;ihit < hits_.size(); ++ihit) {
      auto const& shptr = hits_[ihit];
      auto const& whstate = best[ihit];
      shptr->setState(whstate);
    }
  }

  template<class KTRAJ> void KKStrawHitGroup<KTRAJ>::print(std::ostream& ost, int detail) const {
    unsigned nactive(0), ndrift(0);
    for(auto const& hit : hits_){
      if(hit->active()) ++nactive;
      if(hit->hitState().useDrift()) ++ndrift;
    }
    ost << " KKStrawHitGroup with " << nactive << " active hits with " << ndrift  << " using drift information among " << hits_.size() << " total" << std::endl;
    if(detail > 0){
      for(auto const& hit : hits_) {
        ost << hit->strawId() << std::endl;
      }
    }
  }
}
#endif
