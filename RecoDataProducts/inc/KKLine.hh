#ifndef RecoDataProducts_KKKinematicLine_hh
#define RecoDataProducts_KKKinematicLine_hh
//
// Define a type for storing KinKal Tracks based on KinematicLine
//
#include "KinKal/Fit/Track.hh"
#include "KinKal/Trajectory/KinematicLine.hh"
#include "Offline/GeneralUtilities/inc/OwningPointerCollection.hh"

namespace mu2e {
  using KKLine = KinKal::Track<KinKal::KinematicLine>;
  using KKLineCollection = mu2e::OwningPointerCollection<KKLine>;
}

#endif /* TrkReco_KKLine_hh */
