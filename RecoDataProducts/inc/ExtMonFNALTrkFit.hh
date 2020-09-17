// PatRec results.
//
// Original author Andrei Gaponenko
//
//

#ifndef RecoDataProducts_ExtMonFNALTrkFit_hh
#define RecoDataProducts_ExtMonFNALTrkFit_hh

#include <vector>

#include "canvas/Persistency/Common/Ptr.h"

#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitQuality.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkClusterResiduals.hh"

namespace mu2e {

  class ExtMonFNALRecoCluster;

  //================================================================
  class ExtMonFNALTrkFit {
  public:

    typedef std::vector<art::Ptr<ExtMonFNALRecoCluster> > Clusters;
    typedef std::vector<ExtMonFNALTrkClusterResiduals> Residuals;

    const ExtMonFNALTrkParam& pars() const { return pars_; }
    const ExtMonFNALTrkFitQuality& quality() const { return quality_; }
    const Clusters& clusters() const { return clusters_; }
    const Residuals& residuals() const { return residuals_; }

    ExtMonFNALTrkFit(const ExtMonFNALTrkParam& p,
                     const ExtMonFNALTrkFitQuality& q,
                     const Clusters& c,
                     const Residuals& r)
      : pars_(p), quality_(q), clusters_(c), residuals_(r)
    {}

    // for persistency
    ExtMonFNALTrkFit() {}

  private:
    ExtMonFNALTrkParam pars_;
    ExtMonFNALTrkFitQuality quality_;
    Clusters clusters_;
    Residuals residuals_;
  };

} // namespace mu2e

#endif /* RecoDataProducts_ExtMonFNALTrkFit_hh */
