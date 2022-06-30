
#ifndef RecoDataProducts_KKCentralHelix_hh
#define RecoDataProducts_KKCentralHelix_hh
//
// Define a type for storing KinKal Tracks based on CentralHelix
//
#include "Offline/Mu2eKinKal/inc/KKTrack.hh"
#include "KinKal/Trajectory/CentralHelix.hh"
#include "Offline/GeneralUtilities/inc/OwningPointerCollection.hh"

namespace mu2e {
  using KKCentralHelix = KKTrack<KinKal::CentralHelix>;
  using KKCentralHelixCollection = mu2e::OwningPointerCollection<KKCentralHelix>;
}

#endif /* TrkReco_KKCentralHelix_hh */
