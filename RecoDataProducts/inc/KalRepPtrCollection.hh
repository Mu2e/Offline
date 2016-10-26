#ifndef RecoDataProducts_KalRepPtrCollection_hh
#define RecoDataProducts_KalRepPtrCollection_hh

//
// Define a type for a collection of art::Ptr's to KalRep  objects.
//
// $Id: KalRepPtrCollection.hh,v 1.1 2014/04/18 16:37:04 kutschke Exp $
// $Author: kutschke $
// $Date: 2014/04/18 16:37:04 $
//
// Original author Rob Kutschke
//


// Any class that includes this header needs the following
// using declaration before including this header.  This will
// be needed until the BaBar code is modified.
// using namespace CLHEP:

#include "BTrk/KalmanTrack/KalRep.hh"
#include "canvas/Persistency/Common/Ptr.h"
#include <vector>

namespace mu2e {
  typedef art::Ptr<KalRep>       KalRepPtr;
  typedef std::vector<KalRepPtr> KalRepPtrCollection;
}

#endif /* RecoDataProducts_KalRepPtrCollection_hh */
