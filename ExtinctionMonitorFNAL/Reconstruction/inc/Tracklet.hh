// Andrei Gaponenko, 2012

#ifndef ExtinctionMonitorFNAL_Reconstruction_Tracklet_hh
#define ExtinctionMonitorFNAL_Reconstruction_Tracklet_hh

#include <list>
#include <vector>
#include <ostream>
#include "art/Persistency/Common/Ptr.h"

namespace mu2e {

  class ExtMonFNALRecoCluster;

  namespace ExtMonFNAL {

    //================================================================
    struct Tracklet {

      art::Ptr<ExtMonFNALRecoCluster> firstCluster;
      art::Ptr<ExtMonFNALRecoCluster> lastCluster;

      std::vector<art::Ptr<ExtMonFNALRecoCluster> > middleClusters;

      Tracklet(const art::Ptr<ExtMonFNALRecoCluster>& fc,
               const art::Ptr<ExtMonFNALRecoCluster>& lc)
        : firstCluster(fc)
        , lastCluster(lc)
      {}
    };

    typedef std::list<Tracklet> Tracklets;

    inline std::ostream& operator<<(std::ostream& os, const Tracklet& tl) {
      return os<<"Tracklet(fc="<<tl.firstCluster
               <<", lc="<<tl.lastCluster
               <<", nmiddle="<<tl.middleClusters.size()
               <<" )";
    }

    //================================================================
  } // namespace ExtMonFNAL
} // namespace mu2e

#endif/*ExtinctionMonitorFNAL_Reconstruction_Tracklet_hh*/
