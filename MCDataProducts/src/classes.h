//
// headers needed when genreflex creates root dictionaries
// for objects written to art files
//

#include <vector>
#include <map>

#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Sampled.h"
#include "canvas/Persistency/Provenance/SubRunID.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/RNGsnapshot.h"
#include "cetlib/map_vector.h"

// generation
#include "Offline/MCDataProducts/inc/GenParticle.hh"
#include "Offline/MCDataProducts/inc/GenParticleSPMHistory.hh"
#include "Offline/MCDataProducts/inc/PrimaryParticle.hh"
#include "Offline/MCDataProducts/inc/GenSimParticleLink.hh"
#include "Offline/MCDataProducts/inc/GenEventCount.hh"
#include "Offline/MCDataProducts/inc/StageParticle.hh"

// simulation
#include "Offline/MCDataProducts/inc/StatusG4.hh"
#include "Offline/MCDataProducts/inc/SimParticle.hh"
#include "Offline/MCDataProducts/inc/StepPointMC.hh"
#include "Offline/MCDataProducts/inc/PtrStepPointMCVector.hh"
#include "Offline/MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "Offline/MCDataProducts/inc/SimParticleRemapping.hh"
#include "Offline/MCDataProducts/inc/CosmicLivetime.hh"
#include "Offline/MCDataProducts/inc/SimTimeOffset.hh"
#include "Offline/MCDataProducts/inc/SurfaceStep.hh"

// simulation bookeeping
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfo.hh"
#include "Offline/MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "Offline/MCDataProducts/inc/StepFilterMode.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchTimeMC.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchIntensity.hh"
#include "Offline/MCDataProducts/inc/ProtonBunchIntensitySummary.hh"
#include "Offline/MCDataProducts/inc/EventWeight.hh"

// G4
#include "Offline/MCDataProducts/inc/G4BeamlineInfo.hh"

// MARS
#include "Offline/MCDataProducts/inc/MARSInfo.hh"
#include "Offline/MCDataProducts/inc/GenParticleMARSAssns.hh"
#include "Offline/MCDataProducts/inc/SimParticleMARSAssns.hh"

// calorimeter
#include "Offline/MCDataProducts/inc/CaloMCTruthAssns.hh"
#include "Offline/MCDataProducts/inc/CaloShowerRO.hh"
#include "Offline/MCDataProducts/inc/CaloShowerStep.hh"
#include "Offline/MCDataProducts/inc/CaloShowerSim.hh"
#include "Offline/MCDataProducts/inc/CaloHitMC.hh"
#include "Offline/MCDataProducts/inc/CaloClusterMC.hh"
#include "Offline/MCDataProducts/inc/CaloEDepMC.hh"
// straws
#include "Offline/MCDataProducts/inc/StrawDigiMC.hh"
#include "Offline/MCDataProducts/inc/StrawGasStep.hh"

// tracking
#include "Offline/MCDataProducts/inc/TrackSummaryTruthAssns.hh"
#include "Offline/MCDataProducts/inc/KalSeedMC.hh"

// CRV
#include "Offline/MCDataProducts/inc/CrvStep.hh"
#include "Offline/MCDataProducts/inc/CrvPhotons.hh"
#include "Offline/MCDataProducts/inc/CrvSiPMCharges.hh"
#include "Offline/MCDataProducts/inc/CrvDigiMC.hh"
#include "Offline/MCDataProducts/inc/CrvCoincidenceClusterMC.hh"
#include "Offline/DataProducts/inc/CRSScintillatorBarIndex.hh"
#include "Offline/MCDataProducts/inc/CrvCoincidenceClusterMCAssns.hh"

// ExtMon
#include "Offline/MCDataProducts/inc/ExtMonFNALSimHit.hh"
#include "Offline/MCDataProducts/inc/ExtMonFNALHitTruthAssn.hh"
#include "Offline/MCDataProducts/inc/ExtMonFNALRecoClusterTruthAssn.hh"
#include "Offline/MCDataProducts/inc/ExtMonFNALPatRecTruthAssns.hh"

// Analysis
#include "Offline/MCDataProducts/inc/MCRelationship.hh"

