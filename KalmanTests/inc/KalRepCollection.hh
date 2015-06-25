#ifndef KalmanTests_KalRepCollection_hh
#define KalmanTests_KalRepCollection_hh

//
// Define a type for a collection of KalRep objects.
//
// $Id: KalRepCollection.hh,v 1.1 2012/07/23 17:52:27 brownd Exp $
// $Author: brownd $
// $Date: 2012/07/23 17:52:27 $
//
// Original author Rob Kutschke
//


// Any class that includes this header needs the following
// using declaration before including this header.  This will
// be needed until the BaBar code is modified.
// using namespace CLHEP:

#include "BTrk/KalmanTrack/KalRep.hh"
#include "GeneralUtilities/inc/OwningPointerCollection.hh"

namespace mu2e {
  typedef mu2e::OwningPointerCollection<KalRep> KalRepCollection;
}

#endif /* KalmanTests_KalRepCollection_hh */
