cet_make_library(
    SOURCE
src/CaloClusterMC.cc
src/CaloHitMC.cc
src/CosmicLivetime.cc
src/CrvDigiMC.cc
src/ExtMonFNALSimHit.cc
src/G4BeamlineInfo.cc
src/GenId.cc
src/MARSInfo.cc
src/MCRelationship.cc
src/PhysicalVolumeInfo.cc
src/PhysicalVolumeInfoMultiCollection.cc
src/ProcessCode.cc
src/SimParticle.cc
src/StageParticle.cc
src/StatusG4.cc
src/StepFilterMode.cc
src/StepInstanceName.cc
src/StepPointMC.cc
src/StrawDigiMC.cc
src/TrackSummaryTruthAssns.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
