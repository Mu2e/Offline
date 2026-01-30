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
#include "Offline/RecoDataProducts/inc/RecoProtonBunchIntensity.hh"

// calorimeter
#include "Offline/RecoDataProducts/inc/CaloDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloRecoDigi.hh"
#include "Offline/RecoDataProducts/inc/CaloHit.hh"
#include "Offline/RecoDataProducts/inc/CaloProtoCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloCluster.hh"
#include "Offline/RecoDataProducts/inc/CaloTrigSeed.hh"

// straws
#include "Offline/RecoDataProducts/inc/StrawHit.hh"
#include "Offline/RecoDataProducts/inc/StrawHitFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawFlag.hh"
#include "Offline/RecoDataProducts/inc/StrawDigi.hh"
#include "Offline/RecoDataProducts/inc/StrawDigiFlag.hh"
#include "Offline/RecoDataProducts/inc/ComboHit.hh"

// tracking intermediate products
#include "Offline/RecoDataProducts/inc/HelixHit.hh"
#include "Offline/RecoDataProducts/inc/TimeCluster.hh"
#include "Offline/RecoDataProducts/inc/HelixVal.hh"
#include "Offline/RecoDataProducts/inc/RobustHelix.hh"
#include "Offline/RecoDataProducts/inc/HelixRecoDir.hh"
#include "Offline/RecoDataProducts/inc/HelixSeed.hh"
#include "Offline/RecoDataProducts/inc/BkgCluster.hh"
#include "Offline/RecoDataProducts/inc/BkgClusterHit.hh"
#include "Offline/RecoDataProducts/inc/BkgQual.hh"

// tracking output
#include "Offline/RecoDataProducts/inc/KKLine.hh"

#include "Offline/RecoDataProducts/inc/TrkFitFlag.hh"
#include "Offline/RecoDataProducts/inc/TrkExtTraj.hh"
#include "Offline/RecoDataProducts/inc/KalRepPtrCollection.hh"
#include "Offline/RecoDataProducts/inc/KKLoopHelix.hh"
#include "Offline/RecoDataProducts/inc/KKCentralHelix.hh"
#include "Offline/RecoDataProducts/inc/TrackSummaryRecoMap.hh"
#include "Offline/RecoDataProducts/inc/TrackSummary.hh"
#include "Offline/RecoDataProducts/inc/TrackCaloAssns.hh"
#include "Offline/RecoDataProducts/inc/KalSeed.hh"
#include "Offline/RecoDataProducts/inc/KalIntersection.hh"
#include "Offline/RecoDataProducts/inc/KalSeedAssns.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloHitPID.hh"
#include "Offline/RecoDataProducts/inc/TrkQual.hh"
#include "Offline/RecoDataProducts/inc/RecoQual.hh"
#include "Offline/RecoDataProducts/inc/AlgorithmID.hh"
#include "KinKal/General/ParticleState.hh"
#include "Offline/RecoDataProducts/inc/MVAResult.hh"

// extrapolation and PID
#include "Offline/RecoDataProducts/inc/TrkCaloIntersect.hh"
#include "Offline/RecoDataProducts/inc/TrkCaloMatch.hh"
#include "Offline/RecoDataProducts/inc/PIDProduct.hh"
#include "Offline/RecoDataProducts/inc/TrackClusterMatch.hh"

// CRV
#include "Offline/RecoDataProducts/inc/CrvDigi.hh"
#include "Offline/RecoDataProducts/inc/CrvDAQerror.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulse.hh"
#include "Offline/RecoDataProducts/inc/CrvRecoPulseFlags.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidence.hh"
#include "Offline/RecoDataProducts/inc/CrvCoincidenceCluster.hh"
#include "Offline/RecoDataProducts/inc/CrvStatus.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "artdaq-core-mu2e/Overlays/Decoders/CRVDataDecoder.hh"

// ExtMon
#include "Offline/RecoDataProducts/inc/ExtMonFNALRawHit.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALRawCluster.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALRecoCluster.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALRecoClusterCollection.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALTrkParam.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALTrkFitQuality.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALTrkClusterResiduals.hh"
#include "Offline/RecoDataProducts/inc/ExtMonFNALTrkFit.hh"

// trigger
#include "Offline/RecoDataProducts/inc/TriggerInfo.hh"

// POT / stopped muons monitoring
#include "Offline/RecoDataProducts/inc/IntensityInfoCalo.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTrackerHits.hh"
#include "Offline/RecoDataProducts/inc/IntensityInfoTimeCluster.hh"

// general reco
#include "Offline/RecoDataProducts/inc/RecoCount.hh"

// STM
#include "Offline/RecoDataProducts/inc/STMWaveformDigi.hh"
#include "Offline/RecoDataProducts/inc/STMMWDDigi.hh"
#include "Offline/RecoDataProducts/inc/STMHit.hh"

// MTP
#include "Offline/RecoDataProducts/inc/MTPHit.hh"

#undef ENABLE_MU2E_GENREFLEX_HACKS
