#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"

#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/TrkMaterialTransit.hh"

#include "canvas/Persistency/Common/Assns.h"
namespace mu2e {
  // Assns between a KalSeed and the thing it was created from
  typedef art::Assns<KalSeed,HelixSeed> KalHelixAssns;
  typedef art::Assns<KalSeed,CosmicTrackSeed> KalLineAssns;
  // Assns between a KalSeed and its material transit estimate
  typedef art::Assns<KalSeed,TrkMaterialTransit> KalMaterialTransitAssns;
}
