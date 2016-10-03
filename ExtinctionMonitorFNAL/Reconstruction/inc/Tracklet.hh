// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Reconstruction_Tracklet_hh
#define ExtinctionMonitorFNAL_Reconstruction_Tracklet_hh

#include <list>
#include <vector>
#include <ostream>
#include "canvas/Persistency/Common/Ptr.h"

namespace mu2e {

  class ExtMonFNALRecoCluster;

  namespace ExtMonFNAL {

    //================================================================
    struct Tracklet {

      art::Ptr<ExtMonFNALRecoCluster> firstSeedCluster;
      art::Ptr<ExtMonFNALRecoCluster> secondSeedCluster;

      std::vector<art::Ptr<ExtMonFNALRecoCluster> > addedClusters;

      Tracklet(const art::Ptr<ExtMonFNALRecoCluster>& fc,
               const art::Ptr<ExtMonFNALRecoCluster>& lc)
        : firstSeedCluster(fc)
        , secondSeedCluster(lc)
      {}
    };

    typedef std::list<Tracklet> Tracklets;

    inline std::ostream& operator<<(std::ostream& os, const Tracklet& tl) {
      return os<<"Tracklet(fc="<<tl.firstSeedCluster
               <<", lc="<<tl.secondSeedCluster
               <<", nmiddle="<<tl.addedClusters.size()
               <<" )";
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Reconstruction_Tracklet_hh*/
