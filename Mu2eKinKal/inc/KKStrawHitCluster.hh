#ifndef Mu2eKinKal_KKStrawHitCluster_hh
#define Mu2eKinKal_KKStrawHitCluster_hh
//
//  class describing a cluster of nearby (in time and space) straw hits which need to be updated coherently.
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
  // Functor to define which StrawHits belong in a single cluster
  //
  template <class KTRAJ> class KKStrawHitClusterer {
    public:
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      KKStrawHitClusterer( StrawIdMask clusterLevel, double clusterDt) :
        clusterLevel_(clusterLevel), clusterDt_(clusterDt) {}
      // decide if 2 hits should be clustered
      bool cluster(KKSTRAWHIT const& sh1, KKSTRAWHIT const& sh2) const {
        return clusterLevel_.equal(sh1.strawId(),sh2.strawId()) &&
          fabs(sh1.time()-sh2.time()) < clusterDt_;
      }
      auto clusterLevel() const { return clusterLevel_.level(); }
      double clusterDt() const { return clusterDt_; }
    private:
      StrawIdMask clusterLevel_; // level to cluster straw hits
      double clusterDt_; // max particle time difference between straw hits in a cluster
  };

  template <class KTRAJ> class KKStrawHitCluster : public KinKal::Hit<KTRAJ> {
    public:
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
      using KKSTRAWHITCOL = std::vector<KKSTRAWHITPTR>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using KKSTRAWHITCLUSTERER = KKStrawHitClusterer<KTRAJ>;
      // sort hits by time
      struct StrawHitSort {
        bool operator ()( const KKSTRAWHITPTR& hit1, const KKSTRAWHITPTR& hit2) {
          return hit1->time() < hit2->time(); }
      };
      KKStrawHitCluster() {}
      // create from a single hit
      KKStrawHitCluster(KKSTRAWHITPTR const& hitptr);
      // create from a collection of panel hits
      KKStrawHitCluster(KKSTRAWHITCOL const& hits,KKSTRAWHITCLUSTERER const& clusterer);
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
      ~KKStrawHitCluster(){}
      // KKStrawHitCluster specific interface
      auto const& strawHits() const { return hits_; }
      bool canAddHit(KKSTRAWHITPTR hit,KKSTRAWHITCLUSTERER const& clusterer) const;
      void addHit(KKSTRAWHITPTR hit,KKSTRAWHITCLUSTERER const& clusterer);
    private:
      // references to the individual hits in this hit cluster
      KKSTRAWHITCOL hits_;
  };

  template<class KTRAJ>  KKStrawHitCluster<KTRAJ>::KKStrawHitCluster(KKSTRAWHITPTR const& hitptr) {
    hits_.push_back(hitptr);
  }

  template<class KTRAJ>  KKStrawHitCluster<KTRAJ>::KKStrawHitCluster(KKSTRAWHITCOL const& hits,KKSTRAWHITCLUSTERER const& clusterer) : hits_(hits) {
    // verify clustering
    for (auto ihit=hits_.begin(); ihit != hits_.end(); ++ihit){
      for (auto jhit=ihit+1; jhit != hits_.end(); ++jhit){
        if(!clusterer.cluster(**ihit,**jhit))
          throw cet::exception("RECO")<<"mu2e::KKStrawHitCluster: KKStrawHits don't belong in cluster"<< std::endl;
      }
    }
  }

  template<class KTRAJ> bool KKStrawHitCluster<KTRAJ>::canAddHit(KKSTRAWHITPTR hit, KKSTRAWHITCLUSTERER const& clusterer) const {
    bool add(true);
    for(auto const& myhit : hits_)
      add &= clusterer.cluster(*myhit,*hit);
    return add;
  }

  template<class KTRAJ> void KKStrawHitCluster<KTRAJ>::addHit(KKSTRAWHITPTR hit, KKSTRAWHITCLUSTERER const& clusterer) {
    bool add = canAddHit(hit,clusterer);
    if(!add)throw cet::exception("RECO")<<"mu2e::KKStrawHitCluster: KKStrawHit doesn't belong in cluster"<< std::endl;
    hits_.push_back(hit);
  }

  template<class KTRAJ> double KKStrawHitCluster<KTRAJ>::time() const {
    // return time just past the last hit's time.  This insures hit clusters are updated after individual hits
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

  template<class KTRAJ> void KKStrawHitCluster<KTRAJ>::updateState(KinKal::MetaIterConfig const& miconfig,bool first) {
    if(first){
      // sort the hit ptrs by time
      std::sort(hits_.begin(),hits_.end(),StrawHitSort ());
      // look for an updater; if it's there, update the state
      // Extend this logic if new StrawHitCluster updaters are introduced
      auto cshu = miconfig.findUpdater<CombinatoricStrawHitUpdater>();
      if(cshu != 0){
        cshu->updateHits<KTRAJ>(hits_);
      }
    }
  }

  template<class KTRAJ> void KKStrawHitCluster<KTRAJ>::print(std::ostream& ost, int detail) const {
    unsigned nactive(0), ndrift(0);
    for(auto const& hit : hits_){
      if(hit->active()) ++nactive;
      if(hit->hitState().useDrift()) ++ndrift;
    }
    ost << " KKStrawHitCluster with " << nactive << " active hits with " << ndrift  << " using drift information among " << hits_.size() << " total" << std::endl;
    if(detail > 0){
      for(auto const& hit : hits_) {
        ost << hit->strawId() << std::endl;
      }
    }
  }
}
#endif
