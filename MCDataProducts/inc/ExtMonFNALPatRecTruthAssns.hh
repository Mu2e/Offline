//
// $Id: ExtMonFNALPatRecTruthAssns.hh,v 1.2 2012/11/01 23:38:49 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:38:49 $
//
// Original author Andrei Gaponenko
//

#ifndef RecoDataProducts_ExtMonFNALPatRecTruthAssns_hh
#define RecoDataProducts_ExtMonFNALPatRecTruthAssns_hh

#include "canvas/Persistency/Common/Assns.h"

#include "MCDataProducts/inc/SimParticle.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFit.hh"

namespace mu2e {

  //================================================================
  // Information attached to a particle-trkparam pair in the truth collection

  class ExtMonFNALTrkMatchInfo {
  public:
    explicit ExtMonFNALTrkMatchInfo(unsigned nc, unsigned nt, unsigned np)
      : nCommonClusters_(nc), nTrackClusters_(nt), nParticleClusters_(np)
    {}

    unsigned nCommonClusters()   const { return nCommonClusters_; }
    unsigned nTrackClusters()    const { return nTrackClusters_; }
    unsigned nParticleClusters() const { return nParticleClusters_; }

    // for persistency
    ExtMonFNALTrkMatchInfo()
      : nCommonClusters_(), nTrackClusters_(), nParticleClusters_()
    {}

  private:
    unsigned nCommonClusters_;
    unsigned nTrackClusters_;
    unsigned nParticleClusters_;
  };

  // use as many-to-many Assns
  typedef art::Assns<SimParticle,ExtMonFNALTrkFit,ExtMonFNALTrkMatchInfo> ExtMonFNALPatRecTruthAssns;

} // namespace mu2e

#endif /* RecoDataProducts_ExtMonFNALPatRecTruthAssns_hh */
