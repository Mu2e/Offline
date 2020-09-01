//
// headers needed when genreflex creates root dictionaries
// for objects written to art files
//

#define ENABLE_MU2E_GENREFLEX_HACKS

#include <vector>
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"
#include "RecoDataProducts/inc/CosmicTrack.hh" 
#include "RecoDataProducts/inc/CosmicTrackSeed.hh"

// calorimeter
#include "RecoDataProducts/inc/CaloDigi.hh"
#include "RecoDataProducts/inc/CaloRecoDigiCollection.hh"
#include "RecoDataProducts/inc/CaloHitCollection.hh"
#include "RecoDataProducts/inc/CaloCrystalHit.hh"
#include "RecoDataProducts/inc/CaloProtoClusterCollection.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "RecoDataProducts/inc/CaloTrigSeedCollection.hh"

// straws
#include "RecoDataProducts/inc/StrawHitCollection.hh"
#include "RecoDataProducts/inc/StrawHitFlag.hh"
#include "RecoDataProducts/inc/StrawDigi.hh"
#include "RecoDataProducts/inc/StrawDigiFlag.hh"
#include "RecoDataProducts/inc/ComboHit.hh"

// tracking intermediate products
#include "RecoDataProducts/inc/HelixHit.hh"
#include "RecoDataProducts/inc/TimeCluster.hh"
#include "RecoDataProducts/inc/HelixVal.hh"
#include "RecoDataProducts/inc/RobustHelix.hh"
#include "RecoDataProducts/inc/HelixSeed.hh"
#include "RecoDataProducts/inc/BkgCluster.hh"
#include "RecoDataProducts/inc/BkgClusterHit.hh"
#include "RecoDataProducts/inc/BkgQual.hh"

// tracking output
#include "RecoDataProducts/inc/TrkFitFlag.hh"
#include "RecoDataProducts/inc/TrkExtTrajCollection.hh"
#include "RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "RecoDataProducts/inc/TrackSummaryRecoMap.hh"
#include "RecoDataProducts/inc/TrackSummary.hh"
#include "RecoDataProducts/inc/TrackCaloAssns.hh" 
#include "RecoDataProducts/inc/KalSeed.hh"
#include "RecoDataProducts/inc/TrkCaloHitPID.hh"
#include "RecoDataProducts/inc/TrkQual.hh"
#include "RecoDataProducts/inc/RecoQual.hh"
#include "RecoDataProducts/inc/AlgorithmID.hh"
#include "RecoDataProducts/inc/AlgorithmIDCollection.hh"


// extrapolation and PID
#include "RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "RecoDataProducts/inc/TrkCaloMatchCollection.hh"
#include "RecoDataProducts/inc/AvikPIDProductCollection.hh"
#include "RecoDataProducts/inc/AvikPIDNewProductCollection.hh"
#include "RecoDataProducts/inc/PIDProductCollection.hh"
#include "RecoDataProducts/inc/TrkToCaloExtrapol.hh"
#include "RecoDataProducts/inc/TrackClusterMatch.hh"

// CRV
#include "RecoDataProducts/inc/CrvDigiCollection.hh"
#include "RecoDataProducts/inc/CrvRecoPulseCollection.hh"
#include "RecoDataProducts/inc/CrvCoincidenceCollection.hh"
#include "RecoDataProducts/inc/CrvCoincidenceClusterCollection.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

// ExtMon
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

// trigger
#include "RecoDataProducts/inc/TriggerFlag.hh"
#include "RecoDataProducts/inc/TriggerInfo.hh"

// general reco
#include "RecoDataProducts/inc/RecoCount.hh"

#undef ENABLE_MU2E_GENREFLEX_HACKS
