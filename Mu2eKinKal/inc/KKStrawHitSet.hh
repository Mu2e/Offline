#ifndef Mu2eKinKal_KKStrawHitSet_hh
#define Mu2eKinKal_KKStrawHitSet_hh
//
//  class representing a group of ireferences to nearby (in time and space) straw hits.   This class doesn't constrain
//  the fit, it is used just to coherently update the hits, using correlations to improve the accuracy of ambiguity setting.
//  As such, most of the methods are null-ops; only Update has any purpose
//
#include "KinKal/Detector/Hit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKCombinatoricUpdater.hh"
#include <vector>
namespace mu2e {
  using KinKal::WireHitState;
  template <class KTRAJ> class StrawHitSet : public KinKal::Hit<KTRAJ> {
    public:
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using KKSTRAWHITPTR = shared_ptr<KKSTRAWHIT>;
      // sort hits by time
      struct StrawHitSort {
        bool operator ()( const KKSTRAWHITPTR& hit1, const KKSTRAWHITPTR& hit2) {
          return hit1->time() < hit2->time(); }
      };
      using SHCOLL = std::set<KKSTRAWHITPTR>;
      StrawHitSet() {}
      // create from a collection of panel hits
      StrawHitSet(SHCOLL const& hits) : hits_(hits) {}
      bool addHit(KKSTRAWHITPTR hit); { hits_.insert(hit); }
      Weights weight() const override { return Weights(); }      // Panel hits intrinsically have no weight (null Weight)
      bool active() const override { return false; } // panel hits are never active
      Chisq chisq() const override { return 0.0; } // panel hits don't contribut to chisquared
      Chisq chisq(Parameters const& params) const override { return 0.0; }
      double time() const override;
      void update(PKTRAJ const& pktraj)override {;} // no algebraic update for panel hits
      // update the internals of the hit, specific to this meta-iteraion.  This will affect the next fit iteration
      void update(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      ~KKStrawHitSet(){}
    private:
      // references to the individual hits in this panel hit.  non-const access is needed, to update actual hits
      SHCOLL hits_;
  };

  template<class KTRAJ> bool KKStrawHitSet<KTRAJ>::addHit(KKSTRAWHITPTR hit) {
    hits_.insert(hit);
    return true; // should check that hits are in time FIXME!
  }
  template<class KTRAJ> double KKStrawHitSet<KTRAJ>::time() const {
    // return time just past the last hit's time.  This insures they are updated before the panel
    double maxtime(-std::limits<float>.max());
    unsigned nactive(0);
    for(auto const& hit : hits_){
      if(hit && hit->active())++nactive;
      maxtime = std::max(hit->time(),maxtime);
    }
    static double epsilon(1e-6);
    return maxtime + epsilon;
  }

  template<class KTRAJ> double KKStrawHitSet<KTRAJ>::updateState(PKTRAJ const& pktraj, MetaIterConfig const& config)override {
    // look for an updater; if it's there, update the state
    auto kkphu = miconfig.findUpdater<KKCombinatoricUpdater>();
    if(kkphu != 0)kkphu->update(*this);
  }

  template<class KTRAJ> void KKStrawHitSet<KTRAJ>::print(std::ostream& ost, int detail) const {
    unsigned nactive(0), ndrift(0);
    for(auto const& hit : hits_){
      if(hit->active()) ++nactive;
      if(hit->hitState().useDrift()) ++ndrift;
    }
    ost << " KKStrawHitSet with " << nactive << " active hits with " << ndrift " << using drift information among " << hits_size() << " total" << std::endl;
    if(detail > 0){
      for(auto const& hit : hits_) {
        ost << "hit  << std::endl;
      }
    }
  }
}
#endif
