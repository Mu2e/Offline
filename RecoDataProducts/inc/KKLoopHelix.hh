
#ifndef RecoDataProducts_KKLoopHelix_hh
#define RecoDataProducts_KKLoopHelix_hh
//
// Define a type for storing KinKal Tracks based on LoopHelix
//
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "KinKal/Trajectory/LoopHelix.hh"
#include "Offline/GeneralUtilities/inc/OwningPointerCollection.hh"

namespace mu2e {
  using KKLoopHelix = KKTrack<KinKal::LoopHelix>;
  using KKLoopHelixCollection = mu2e::OwningPointerCollection<KKLoopHelix>;
}

#endif /* TrkReco_KKLoopHelix_hh */
