//
// headers needed when genreflex creates root dictionaries
// for objects written to art files
//

#define ENABLE_MU2E_GENREFLEX_HACKS

#include <vector>
#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Assns.h"
#include "Offline/RecoDataProducts/inc/CosmicTrack.hh" 
#include "Offline/RecoDataProducts/inc/CosmicTrackSeed.hh"

// beam
#include "Offline/RecoDataProducts/inc/ProtonBunchTime.hh"

// calorimeter
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloProtoCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloTrigSeed.hh"

// straws
#include "Offline/RecoDataProducts/inc/StrawHitCollection.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawDigiFlag.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

// tracking intermediate products
#include "Offline/RecoDataProducts/inc/HelixHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/HelixVal.hh"
#include "Offline/RecoDataProducts/inc/RobustHelix.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/BkgCluster.hh"
#include "Offline/RecoDataProducts/inc/BkgClusterHit.hh"
#include "Offline/RecoDataProducts/inc/BkgQual.hh"

// tracking output
#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "Offline/RecoDataProducts/inc/TrkExtTrajCollection.hh"
#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "Offline/RecoDataProducts/inc/KKLoopHelix.hh"
#include "Offline/RecoDataProducts/inc/TrackSummaryRecoMap.hh"
#include "Offline/RecoDataProducts/inc/TrackSummary.hh"
#include "Offline/RecoDataProducts/inc/TrackCaloAssns.hh" 
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloHitPID.hh"
#include "Offline/RecoDataProducts/inc/TrkQual.hh"
#include "Offline/RecoDataProducts/inc/RecoQual.hh"
#include "Offline/RecoDataProducts/inc/AlgorithmID.hh"
#include "Offline/RecoDataProducts/inc/AlgorithmIDCollection.hh"
#include "KinKal/General/ParticleState.hh"

// extrapolation and PID
#include "Offline/RecoDataProducts/inc/TrkCaloIntersectCollection.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloMatchCollection.hh"
#include "Offline/RecoDataProducts/inc/AvikPIDProductCollection.hh"
#include "Offline/RecoDataProducts/inc/AvikPIDNewProductCollection.hh"
#include "Offline/RecoDataProducts/inc/PIDProductCollection.hh"
#include "Offline/RecoDataProducts/inc/TrkToCaloExtrapol.hh"
#include "Offline/RecoDataProducts/inc/TrackClusterMatch.hh"

// CRV
#include "Offline/RecoDataProducts/inc/CrvDigiCollection.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulseFlags.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCollection.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceClusterCollection.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"

// ExtMon
#include "Offline/RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALRawHitCollection.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALRawCluster.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALRawClusterCollection.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALTrkFitQuality.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALTrkClusterResiduals.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALTrkFit.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALTrkFitCollection.hh"

// trigger
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"

// POT / stopped muons monitoring bvitaly May 2021
#include "Offline/RecoDataProducts/inc/IntensityInfo.hh"

// general reco
#include "Offline/RecoDataProducts/inc/RecoCount.hh"

#undef ENABLE_MU2E_GENREFLEX_HACKS
