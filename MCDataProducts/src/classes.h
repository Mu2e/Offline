//
// headers needed when genreflex creates root dictionaries
// for objects written to art files
//

#include <vector>
#include <map>

#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/RNGsnapshot.h"
#include "cetlib/map_vector.h"

// generation
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/GenParticleCollections.hh"
#include "MCDataProducts/inc/GenParticleSPMHistory.hh"
#include "MCDataProducts/inc/PrimaryParticle.hh"
#include "MCDataProducts/inc/GenSimParticleLink.hh"
#include "MCDataProducts/inc/GenEventCount.hh"
#include "MCDataProducts/inc/FixedTimeMap.hh"
#include "MCDataProducts/inc/StageParticle.hh"

// simulation
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"
#include "MCDataProducts/inc/CosmicLivetime.hh"
#include "MCDataProducts/inc/SimTimeOffset.hh"

// simulation bookeeping
#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "MCDataProducts/inc/StepFilterMode.hh"
#include "MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "MCDataProducts/inc/ProtonBunchIntensitySummary.hh"
#include "MCDataProducts/inc/EventWeight.hh"

// G4
#include "MCDataProducts/inc/G4BeamlineInfo.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"

// MARS
#include "MCDataProducts/inc/MARSInfo.hh"
#include "MCDataProducts/inc/MARSInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleMARSAssns.hh"
#include "MCDataProducts/inc/SimParticleMARSAssns.hh"

// calorimeter
#include "MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "MCDataProducts/inc/CaloShowerRO.hh"
#include "MCDataProducts/inc/CaloShowerStep.hh"
#include "MCDataProducts/inc/CaloShowerSim.hh"
#include "MCDataProducts/inc/CaloHitMC.hh"
#include "MCDataProducts/inc/CaloClusterMC.hh"
#include "MCDataProducts/inc/CaloEDepMC.hh"
// straws
#include "MCDataProducts/inc/StrawDigiMC.hh"
#include "MCDataProducts/inc/StrawGasStep.hh"

// tracking 
#include "MCDataProducts/inc/TrackSummaryTruthAssns.hh"
#include "MCDataProducts/inc/KalSeedMC.hh"

// CRV
#include "MCDataProducts/inc/CrvStep.hh"
#include "MCDataProducts/inc/CrvPhotons.hh"
#include "MCDataProducts/inc/CrvSiPMCharges.hh"
#include "MCDataProducts/inc/CrvDigiMC.hh"
#include "MCDataProducts/inc/CrvCoincidenceClusterMCCollection.hh"
#include "DataProducts/inc/CRSScintillatorBarIndex.hh"

// ExtMon
#include "MCDataProducts/inc/ExtMonFNALSimHit.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALHitTruthAssn.hh"
#include "MCDataProducts/inc/ExtMonFNALRecoClusterTruthAssn.hh"
#include "MCDataProducts/inc/ExtMonFNALPatRecTruthAssns.hh"

// Analysis
#include "MCDataProducts/inc/MCRelationship.hh"

