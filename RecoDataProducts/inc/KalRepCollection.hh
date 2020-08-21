#ifndef TrkReco_KalRepCollection_hh
#define TrkReco_KalRepCollection_hh
//
// Define a type for a collection of KalRep objects.
//
//
// Original author Rob Kutschke
//
#include "BTrk/KalmanTrack/KalRep.hh"
#include "GeneralUtilities/inc/OwningPointerCollection.hh"

namespace mu2e {
  typedef mu2e::OwningPointerCollection<KalRep> KalRepCollection;
}

#endif /* TrkReco_KalRepCollection_hh */
