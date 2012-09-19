//
// $Id: ExtMonFNALPatRecTrackAssns.hh,v 1.1 2012/09/19 03:54:19 gandr Exp $
// $Author: gandr $
// $Date: 2012/09/19 03:54:19 $
//
// Original author Andrei Gaponenko
//

#ifndef RecoDataProducts_ExtMonFNALPatRecTrackAssns_hh
#define RecoDataProducts_ExtMonFNALPatRecTrackAssns_hh


#include "art/Persistency/Common/Assns.h"

#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"

namespace mu2e {

  // use as a many-to-many Assns
  typedef art::Assns<ExtMonFNALRecoCluster,ExtMonFNALTrkParam> ExtMonFNALPatRecTrackAssns;

} // namespace mu2e

#endif /* RecoDataProducts_ExtMonFNALPatRecTrackAssns_hh */
