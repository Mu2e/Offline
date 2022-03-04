#ifndef Mu2eKinKal_KKPanelHit_hh
#define Mu2eKinKal_KKPanelHit_hh
//
//  class representing a set of straw hits from a single particle in a single panel.  This class doesn't constrain
//  the fit, it is used just to coherently update the hits, using correlations to improve the accuracy of ambiguity setting.
//  As such, most of the methods are null-ops; only Update has any purpose
//
#include "KinKal/Detector/Hit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include "Offline/Mu2eKinKal/inc/KKPanelHitUpdater.hh"
#include <vector>
namespace mu2e {
  using KinKal::WireHitState;
  template <class KTRAJ> class PanelHit : public KinKal::Hit<KTRAJ> {
    public:
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using SHCOLL = std::vector<shared_ptr<KKSTRAWHIT>>;
      // create from a collection of panel hits
      PanelHit(SHCOLL const& hits) : hits_(hits) {}
      Weights weight() const override { return Weights(); }      // Panel hits intrinsically have no weight
      bool active() const override { return false; } // panel hits are never active
      Chisq chisq() const override { return 0.0; } // panel hits don't contribut to chisquared
      Chisq chisq(Parameters const& params) const override { return 0.0; }
      double time() const override;
      void update(PKTRAJ const& pktraj)override {;} // no algebraic update for panel hits
      // update the internals of the hit, specific to this meta-iteraion.  This function is the reason for this class
      void update(PKTRAJ const& pktraj, MetaIterConfig const& config) override;
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      ~KKPanelHit(){}
    private:
      // references to the individual hits in this panel hit.  non-const access is needed, to update actual hits
      SHCOLL hits_;
  };

  template<class KTRAJ> double KKPanelHit<KTRAJ>::time() const {
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

  template<class KTRAJ> double KKPanelHit<KTRAJ>::updateState(PKTRAJ const& pktraj, MetaIterConfig const& config)override {
  // look for an updater; if it's there, update the state
    auto kkphu = miconfig.findUpdater<KKPanelHitUpdater>();
    if(kkphu != 0)kkphu->update(*this);
  }

  template<class KTRAJ> void KKPanelHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    unsigned nactive(0), ndrift(0);
    for(auto const& hit : hits_){
      if(hit->active()) ++nactive;
      if(hit->hitState().useDrift()) ++ndrift;
    }
    ost << " KKPanelHit with " << nactive << " active hits with " << ndrift " << using drift information among " << hits_size() << " total" << std::endl;
    if(detail > 0){
      for(auto const& hit : hits_) {
          ost << "hit  << std::endl;
      }
    }
  }
}
#endif
