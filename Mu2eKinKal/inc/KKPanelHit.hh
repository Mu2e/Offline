#ifndef Mu2eKinKal_KKPanelHit_hh
#define Mu2eKinKal_KKPanelHit_hh
//
//  class representing a set of straw hits from a single particle in a single panel.  This class doesn't constrain
//  the fit, it is used just to coherently upcate the hits, using correlations to improve the accuracy of ambiguity setting.
//
#include "KinKal/Detector/Hit.hh"
#include "Offline/Mu2eKinKal/inc/KKStrawHit.hh"
#include <array>
namespace mu2e {
  // struct for updating straw hits;
  struct KKPanelHitUpdater {
    double mindoca_; // minimum DOCA value to set an ambiguity
    double maxdoca_; // maximum DOCA to still use a hit
    bool nulltime_; // constrain time when hit has null ambiguity
    double rcell_; // straw radius
    KKStrawHitUpdater(double mindoca,double maxdoca, bool nulltime) : mindoca_(mindoca), maxdoca_(maxdoca), nulltime_(nulltime) {}
  };

  using KinKal::WireHitState;
  template <class KTRAJ> class PanelHit : public KinKal::Hit<KTRAJ> {
    public:
      static constexpr MAXNHIT=4;
      using KKSTRAWHIT = KKStrawHit<KTRAJ>;
      using SHCOLL = std::array<shared_ptr<KKSTRAWHIT>,MAXNHIT>;
      using PHState = std::array<WireHitState,MAXNHIT>;
// the constraint this hit implies WRT the current reference, expressed as a weight
      Weights weight() const override;
      // hits are active if any component is active
      bool active() const override { return false; } // this hit never contributes to the fit
      Chisq chisq() const override { return 0.0; }
      Chisq chisq(Parameters const& params) const override { return 0.0; }
      double time() const override;
      void update(PKTRAJ const& pktraj)override {;} // no POCA update for this class
      // update the internals of the hit, specific to this meta-iteraion.  This function is the reason for this class
      void updateState(PKTRAJ const& pktraj, MetaIterConfig const& config)override;
      PanelHit(SHCOLL const& hits) : hits_(hits) {}
      void print(std::ostream& ost=std::cout,int detail=0) const override;
      ~KKPanelHit(){}
    private:
      // references to the individual hits in this panel.  note non-const access is needed
      SHCOLL hits_;
      PHState pstate_; // aggregate state of this panel's hits
  };

  template<class KTRAJ> double KKPanelHit<KTRAJ>::time() const {
    // return average of hits
    double time(0.0);
    unsigned nactive(0);
    for(auto const& hit : hits_){
      if(hit && hit->active()){
  ++nactive;
  time += hit->time();
      }
    }
    if(nactive > 0)
      time /= static_cast<float>(nactive);
    else  // no active hits: just take the 1st one
      time = hits_.front()->time();
    return time;
  }

  template<class KTRAJ> double KKPanelHit<KTRAJ>::updateState(PKTRAJ const& pktraj, MetaIterConfig const& config)override {
// find the specific configuration for this meta-iteration

  }

  template<class KTRAJ> void KKPanelHit<KTRAJ>::print(std::ostream& ost, int detail) const {
    unsigned nhits(0), nactive(0);
    for(auto const& hit : hits_){
      if(hit){
  ++nhits;
  if(hit->active()){
    ++nactive;
  }
      }
    }
    ost << " KKPanelHit with " << nactive << " active hits among " << nhits << " total" << std::endl;
    if(detail > 0){
      for(size_t ihit=0; ihit < MAXNHIT; ++ihit){
  if(hits_[ihit]){
    ost << pstate_[ihit] << std::endl;
  }
      }
    }
  }
}
#endif
