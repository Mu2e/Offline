//
// Original author Rob Kutschke
//

#define ENABLE_MU2E_GENREFLEX_HACKS

#include <vector>

#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"

#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHitCollection.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"
#include "RecoDataProducts/inc/CaloClusterCollection.hh"
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StereoHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitPositionCollection.hh"
#include "RecoDataProducts/inc/HelixHit.hh"
#include "RecoDataProducts/inc/StrawHitFlagCollection.hh"
#include "RecoDataProducts/inc/StrawClusterCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulsesCollection.hh"
#include "RecoDataProducts/inc/CrvCoincidenceCheckResult.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/TrackSummaryRecoMap.hh"
#include "RecoDataProducts/inc/ExtMonUCITofHitCollection.hh"
#include "RecoDataProducts/inc/SubEventCollection.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/HelixVal.hh"
#include "RecoDataProducts/inc/RobustHelix.hh"
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "RecoDataProducts/inc/ZRotStrawHitMapCollection.hh"
#include "RecoDataProducts/inc/TrackerHitByID.hh"
#include "RecoDataProducts/inc/TrackSeedCollection.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/TrkExtTrajCollection.hh"
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrkCaloMatchCollection.hh"
#include "RecoDataProducts/inc/AvikPIDProductCollection.hh"
#include "RecoDataProducts/inc/AvikPIDNewProductCollection.hh"
#include "RecoDataProducts/inc/PIDProductCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRawClusterCollection.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"

#include "RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitQuality.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkClusterResiduals.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFit.hh"
#include "RecoDataProducts/inc/ExtMonFNALTrkFitCollection.hh"

#include "RecoDataProducts/inc/StrawDigiCollection.hh"
#include "RecoDataProducts/inc/TrackSummary.hh"

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "RecoDataProducts/inc/TrackCaloAssns.hh" 
// BTrk Kalman fit persistence
#include "RecoDataProducts/inc/KalSeed.hh"

#include "RecoDataProducts/inc/TriggerFlag.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"

#include "RecoDataProducts/inc/TrkQual.hh"

#undef ENABLE_MU2E_GENREFLEX_HACKS
