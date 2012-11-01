// PatRec results.
//
// Original author Andrei Gaponenko
//
//
// $Id: ExtMonFNALTrkFit.hh,v 1.1 2012/11/01 23:38:21 gandr Exp $
// $Author: gandr $
// $Date: 2012/11/01 23:38:21 $

#ifndef RecoDataProducts_ExtMonFNALTrkFit_hh
#define RecoDataProducts_ExtMonFNALTrkFit_hh

#include <vector>

#include "art/Persistency/Common/Ptr.h"

#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitQuality.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkClusterResiduals.hh"

namespace mu2e {

  class ExtMonFNALRecoCluster;

  //================================================================
  struct ExtMonFNALTrkFit {
    ExtMonFNALTrkParam pars;
    ExtMonFNALTrkFitQuality quality;
    std::vector<art::Ptr<ExtMonFNALRecoCluster> > clusters;
    std::vector<ExtMonFNALTrkClusterResiduals> residuals;
  };

} // namespace mu2e

#endif /* RecoDataProducts_ExtMonFNALTrkFit_hh */
