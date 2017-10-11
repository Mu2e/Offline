//
// Original author Rob Kutschke
//

//
// For every type that is described in classes_def.xml, the header must
// be included in this file. It is OK if the header comes in as the
// result of including another header.
//
// Starting with ROOT6, the template instantiations should no longer be
// present.
//
#include <vector>
#include <map>

#include "canvas/Persistency/Common/Wrapper.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/Assns.h"
#include "canvas/Persistency/Common/RNGsnapshot.h"
#include "cetlib/map_vector.h"

#include "MCDataProducts/inc/CaloHitMCTruthAssn.hh"
#include "MCDataProducts/inc/CaloClusterMCTruthAssn.hh"
#include "MCDataProducts/inc/CaloShowerStepROCollection.hh"
#include "MCDataProducts/inc/CaloShowerStepCollection.hh"
#include "MCDataProducts/inc/CaloShowerSimCollection.hh"
#include "RecoDataProducts/inc/CaloDigiCollection.hh"
#include "MCDataProducts/inc/CaloHitMCTruthCollection.hh"
#include "MCDataProducts/inc/CaloDigiMCCollection.hh"
#include "MCDataProducts/inc/CaloHitSimPartMCCollection.hh"
#include "MCDataProducts/inc/CaloClusterContentMC.hh"
#include "RecoDataProducts/inc/CaloCluster.hh"
#include "MCDataProducts/inc/CrvPhotonArrivalsCollection.hh"
#include "MCDataProducts/inc/CrvSiPMResponsesCollection.hh"
#include "MCDataProducts/inc/CrvDigiMCCollection.hh"
#include "MCDataProducts/inc/ExtMonUCITofHitMCTruthCollection.hh"
#include "MCDataProducts/inc/G4BeamlineInfo.hh"
#include "MCDataProducts/inc/G4BeamlineInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfo.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoCollection.hh"
#include "MCDataProducts/inc/PhysicalVolumeInfoMultiCollection.hh"
#include "MCDataProducts/inc/MCTrajectoryCollection.hh"
#include "MCDataProducts/inc/PointTrajectoryCollection.hh"
#include "MCDataProducts/inc/SimParticleCollection.hh"
#include "MCDataProducts/inc/SimParticlePtrCollection.hh"
#include "MCDataProducts/inc/StatusG4.hh"
#include "MCDataProducts/inc/StepPointMCCollection.hh"
#include "MCDataProducts/inc/StrawHitMCTruthCollection.hh"
#include "MCDataProducts/inc/StrawDigiMCCollection.hh"
#include "MCDataProducts/inc/PtrStepPointMCVectorCollection.hh"
#include "MCDataProducts/inc/MixingSummary.hh"
#include "MCDataProducts/inc/VisibleGenElTrackCollection.hh"
#include "MCDataProducts/inc/GenParticleSPMHistory.hh"
#include "MCDataProducts/inc/GenSimParticleLink.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHit.hh"
#include "MCDataProducts/inc/ExtMonFNALSimHitCollection.hh"
#include "MCDataProducts/inc/ExtMonFNALHitTruthAssn.hh"
#include "MCDataProducts/inc/ExtMonFNALRecoClusterTruthAssn.hh"
#include "MCDataProducts/inc/ExtMonFNALPatRecTruthAssns.hh"
#include "MCDataProducts/inc/MARSInfo.hh"
#include "MCDataProducts/inc/MARSInfoCollection.hh"
#include "MCDataProducts/inc/GenParticleMARSAssns.hh"
#include "MCDataProducts/inc/SimParticleMARSAssns.hh"
#include "MCDataProducts/inc/SimParticleTimeMap.hh"
#include "MCDataProducts/inc/SimParticleRemapping.hh"

#include "MCDataProducts/inc/StepFilterMode.hh"
#include "MCDataProducts/inc/GenEventCount.hh"
#include "MCDataProducts/inc/EventWeight.hh"

#include "MCDataProducts/inc/TrackSummaryTruthAssns.hh"
#include "MCDataProducts/inc/SimParticleTrackerHitSummaryCollection.hh"
#include "MCDataProducts/inc/SimParticleTrackerHitSummaryAssns.hh"
#include "MCDataProducts/inc/ProtonBunchIntensity.hh"

#include "DataProducts/inc/CRSScintillatorBarIndex.hh"
