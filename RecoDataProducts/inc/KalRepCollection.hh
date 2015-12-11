#ifndef TrkReco_KalRepCollection_hh
#define TrkReco_KalRepCollection_hh
//
// Define a type for a collection of KalRep objects.
//
// $Id: KalRepCollection.hh,v 1.1 2012/07/23 17:52:27 brownd Exp $
// $Author: brownd $
// $Date: 2012/07/23 17:52:27 $
//
// Original author Rob Kutschke
//
#include "BTrk/KalmanTrack/KalRep.hh"
#include "GeneralUtilities/inc/OwningPointerCollection.hh"

namespace mu2e {
  typedef mu2e::OwningPointerCollection<KalRep> KalRepCollection;
}

#endif /* TrkReco_KalRepCollection_hh */
