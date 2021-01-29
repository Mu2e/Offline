
#ifndef RecoDataProducts_KKLoopHelix_hh
#define RecoDataProducts_KKLoopHelix_hh
//
// Define a type for storing KinKal Tracks based on LoopHelix
//
#include "KinKal/Fit/Track.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "GeneralUtilities/inc/OwningPointerCollection.hh"

namespace mu2e {
  using KKLoopHelix = KinKal::Track<KinKal::LoopHelix>;
  using KKLoopHelixCollection = mu2e::OwningPointerCollection<KKLoopHelix>;
}

#endif /* TrkReco_KKLoopHelix_hh */
