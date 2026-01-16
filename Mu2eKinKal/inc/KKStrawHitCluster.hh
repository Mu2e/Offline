#ifndef Mu2eKinKal_KKStrawHitCluster_hh
#define Mu2eKinKal_KKStrawHitCluster_hh
//
//  class describing a cluster of nearby (in time and space) straw hits which need to be updated coherently.
//  Even though it is formally a 'Hit', its only interaction with the fit is in iteration updating.
//
#include "KinKal/Detector/Hit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawXing.hh"
#include "Offline/Mu2eKinKal/inc/WHSIterator.hh"
#include "Offline/Mu2eKinKal/inc/Chi2SHU_updateCluster.hh"
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
      KKStrawHitClusterer( StrawIdMask clusterLevel, size_t maxdstraw, double clusterDt) :
        clusterLevel_(clusterLevel), maxdstraw_(maxdstraw), clusterDt_(clusterDt) {}
      // decide if 2 hits should be clustered
      bool cluster(KKSTRAWHIT const& sh1, KKSTRAWHIT const& sh2) const {
        return clusterLevel_.equal(sh1.strawId(),sh2.strawId()) &&
          abs(sh1.strawId().straw()-sh2.strawId().straw())<= maxdstraw_ &&
          fabs(sh1.time()-sh2.time()) < clusterDt_;
      }
      auto clusterLevel() const { return clusterLevel_.level(); }
      double clusterDt() const { return clusterDt_; }
    private:
      StrawIdMask clusterLevel_; // level to cluster straw hits
      size_t maxdstraw_; // max straw index difference
      double clusterDt_; // max particle time difference between straw hits in a cluster
  };

  template <class KTRAJ> class KKStrawHitCluster : public KinKal::Hit<KTRAJ> {
    public:
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using KKSTRAWHITPTR = std::shared_ptr<KKSTRAWHIT>;
      using KKSTRAWHITCOL = std::vector<KKSTRAWHITPTR>;
      using KKSTRAWXING = KKStrawXing<KTRAJ>;
      using KKSTRAWXINGPTR = std::shared_ptr<KKSTRAWXING>;
      using KKSTRAWXINGCOL = std::vector<KKSTRAWXINGPTR>;
      using KTRAJPTR = std::shared_ptr<KTRAJ>;
      using KKSTRAWHITCLUSTERER = KKStrawHitClusterer<KTRAJ>;
      using PTRAJ = KinKal::ParticleTrajectory<KTRAJ>;
      KKStrawHitCluster() {}
      // create from a single hit
      KKStrawHitCluster(KKSTRAWHITPTR const& hitptr);
      // create from a collection of panel hits
      KKStrawHitCluster(KKSTRAWHITCOL const& hits,KKSTRAWHITCLUSTERER const& clusterer);
      // clone op for reinstantiation
      KKStrawHitCluster(KKStrawHitCluster<KTRAJ> const& rhs){
        /**/
      };
      std::shared_ptr< KinKal::Hit<KTRAJ> > clone(CloneContext& context) const override{
        auto rv = std::make_shared< KKStrawHitCluster<KTRAJ> >(*this);
        for (const auto& ptr: this->strawHits()){
          auto hit = context.get(ptr);
          rv->push_back(hit);
        }
        for (const auto& ptr: this->strawXings()){
          auto xng = context.get(ptr);
          rv->push_back(xng);
        }
        return rv;
      };
      //Hit interface
      bool active() const override { return false; } // panel hits are never active
      KinKal::Chisq chisq(KinKal::Parameters const& params) const override { return KinKal::Chisq(); }
      unsigned nDOF() const override { return 0; }
      KinKal::Weights const& weight() const override { return (*hits_.begin())->weight(); }
      double time() const override;
      void updateReference(PTRAJ const& ptraj) override {} // nothing to do here, ref comes from individual hits
      KTRAJPTR refTrajPtr() const override { return (*hits_.begin())->refTrajPtr(); }
      // update the internals of the hit, specific to this meta-iteraion.  This will affect the next fit iteration
      void updateState(KinKal::MetaIterConfig const& config,bool first) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      ~KKStrawHitCluster(){}
      // KKStrawHitCluster specific interface
      KinKal::TimeRange timeRange() const;
      auto const& strawHits() const { return hits_; }
      auto const& strawXings() const { return xings_; }
      bool canAddHit(KKSTRAWHITPTR hit,KKSTRAWHITCLUSTERER const& clusterer) const;
      void addHit(KKSTRAWHITPTR hit,KKSTRAWHITCLUSTERER const& clusterer);
      void addXing(KKSTRAWXINGPTR xing);

    protected:
      void push_back(KKSTRAWHITPTR hit){ hits_.push_back(hit); }
      void push_back(KKSTRAWXINGPTR xng){ xings_.push_back(xng); };

    private:
      // references to the individual hits and xings in this hit cluster
      KKSTRAWHITCOL hits_;
      KKSTRAWXINGCOL xings_;
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

  template<class KTRAJ> void KKStrawHitCluster<KTRAJ>::addXing(KKSTRAWXINGPTR sxing) {
// make sure this Xing isn't already here
    bool exists(false);
    for(auto const& xing : xings_) {
      if(sxing->strawId() == xing->strawId()){
        exists = true;
        break;
      }
    }
    if(!exists)xings_.push_back(sxing);
  }

  template<class KTRAJ> double KKStrawHitCluster<KTRAJ>::time() const {
    // return time just before the first hit's time.  This insures hit clusters are updated before individual hits
    // This insures the weights subtracted correspond to the reference fit, and that any changes made to the
    // hits in the cluster get propagated to the residuals and weights before the next fit
    double mintime(std::numeric_limits<float>::max());
    for(auto const& hit : hits_){
      mintime = std::min(hit->time(),mintime);
    }
    static double epsilon(1e-6); // small time
    return mintime - epsilon;
  }

  template<class KTRAJ> KinKal::TimeRange KKStrawHitCluster<KTRAJ>::timeRange() const {
    double tmin(0.0), tmax(0.0);
   if(hits_.size() >0){
      tmin = tmax = hits_.front()->time();
      for (auto const& hit: hits_) {
        tmin = std::min(hit->time(),tmin);
        tmax = std::max(hit->time(),tmax);
      }
    }
    return KinKal::TimeRange(tmin,tmax);
  }

  template<class KTRAJ> void KKStrawHitCluster<KTRAJ>::updateState(KinKal::MetaIterConfig const& miconfig,bool first) {
    if(first){
      // look for an updater; if it's there, update the state
      // Extend this logic if new StrawHitCluster updaters are introduced
      auto cshu = miconfig.findUpdater<Chi2SHU>();
      if(cshu != 0){
        cshu->updateCluster<KTRAJ>(*this,miconfig);
      }
    }
    // on subsequent iterations the cluster is left unchanged (but the residuals still are updated)
  }

  template<class KTRAJ> void KKStrawHitCluster<KTRAJ>::print(std::ostream& ost, int detail) const {
    unsigned nactive(0), ndrift(0);
    for(auto const& hit : hits_){
      if(hit->active()) ++nactive;
      if(hit->hitState().driftConstraint()) ++ndrift;
    }
    ost << " KKStrawHitCluster with " << nactive << " active hits with " << ndrift  << " using drift information among " << hits_.size() << " total" << std::endl;
    if(detail > 0){
      ost << "Hits ";
      for(auto const& hit : hits_) ost << hit->strawId() << " time " << hit->time();
      ost << " Xings ";
      for(auto const& xing : xings_) ost << xing->strawId() << " time " << xing->time() << std::endl;
    }
  }
}
#endif
