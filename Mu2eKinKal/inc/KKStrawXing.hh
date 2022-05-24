#ifndef Mu2eKinKal_KKStrawXing_hh
#define Mu2eKinKal_KKStrawXing_hh
//
//  Subclass of StrawXing including Mu2e-specific payload
//
#include "KinKal/Detector/StrawXing.hh"
#include "KinKal/Trajectory/PiecewiseClosestApproach.hh"
#include "Offline/DataProducts/inc/StrawId.hh"
namespace mu2e {
  template <class KTRAJ> class KKStrawXing : public KinKal::StrawXing<KTRAJ> {
    public:
      using PTCA = KinKal::PiecewiseClosestApproach<KTRAJ,Line>;
      KKStrawXing(PTCA const& ptca, KinKal::StrawMaterial const& smat, StrawId sid) : KinKal::StrawXing<KTRAJ>(ptca,smat), sid_(sid) {}
      StrawId strawId() const { return sid_; }
    private:
      StrawId sid_;
  };
}
#endif
