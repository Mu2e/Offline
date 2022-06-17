//
// KinKal fit module for LoopHelix
//
// Original author D. Brown (LBNL) 11/18/20
//
#include "KinKal/Trajectory/LoopHelix.hh"
using KTRAJ= KinKal::LoopHelix;
#include "Offline/Mu2eKinKal/inc/HelixFit_module.hh"
namespace mu2e {
  class LoopHelixFit : public HelixFit {
    public:
      explicit LoopHelixFit(const GlobalSettings& settings) :
        HelixFit(settings,TrkFitFlag::KKLoopHelix) {}
      virtual ~LoopHelixFit() {}
  };
}
DEFINE_ART_MODULE(mu2e::LoopHelixFit);
