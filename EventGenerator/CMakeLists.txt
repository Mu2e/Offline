cet_make_library(
    SOURCE
src/CosmicCRY.cc
src/CosmicDYB.cc
src/CosmicFromTH2.cc
src/DIOGenerator_tool.cc
src/DYBGenerator.cc
src/ExtMonFNALGunImpl.cc
src/FromG4BLFile.cc
src/MuCap1809keVGammaGenerator_tool.cc
src/MuCapDeuteronGenerator_tool.cc
src/MuCapNeutronGenerator_tool.cc
src/MuCapPhotonGenerator_tool.cc
src/MuCapProtonGenerator_tool.cc
src/MuplusMichelGenerator_tool.cc
src/ParticleGun.cc
src/ParticleGunImpl.cc
src/PiCaptureEffects.cc
src/PrimaryProtonGunImpl.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(AntiProtonGun art::module REG_SOURCE src/AntiProtonGun_module.cc)

cet_build_plugin(CaloCalibGun art::module REG_SOURCE src/CaloCalibGun_module.cc)

cet_build_plugin(CaloTBGun art::module REG_SOURCE src/CaloTBGun_module.cc)

cet_build_plugin(CeEndpoint art::module REG_SOURCE src/CeEndpoint_module.cc)

cet_build_plugin(CORSIKAEventGenerator art::module REG_SOURCE src/CORSIKAEventGenerator_module.cc)

cet_build_plugin(CRYEventGenerator art::module REG_SOURCE src/CRYEventGenerator_module.cc)

cet_build_plugin(CryResampler art::module REG_SOURCE src/CryResampler_module.cc)

cet_build_plugin(DecayInOrbitWeight art::module REG_SOURCE src/DecayInOrbitWeight_module.cc)

cet_build_plugin(EplusFromStoppedPion art::module REG_SOURCE src/EplusFromStoppedPion_module.cc)

cet_build_plugin(EventGenerator art::module REG_SOURCE src/EventGenerator_module.cc)

cet_build_plugin(ExtMonFNALBoxGenerator art::module REG_SOURCE src/ExtMonFNALBoxGenerator_module.cc)

cet_build_plugin(ExtMonFNALGun art::module REG_SOURCE src/ExtMonFNALGun_module.cc)

cet_build_plugin(ExtMonFNALRoomGenerator art::module REG_SOURCE src/ExtMonFNALRoomGenerator_module.cc)

cet_build_plugin(FlatMuonDaughterGenerator art::module REG_SOURCE src/FlatMuonDaughterGenerator_module.cc)

cet_build_plugin(FromAsciiMomentumAndPosition art::module REG_SOURCE src/FromAsciiMomentumAndPosition_module.cc)

cet_build_plugin(FromSimParticleCompact art::module REG_SOURCE src/FromSimParticleCompact_module.cc)

cet_build_plugin(FromStepPointMCs art::module REG_SOURCE src/FromStepPointMCs_module.cc)

cet_build_plugin(FromStepPointMCsRotateTarget art::module REG_SOURCE src/FromStepPointMCsRotateTarget_module.cc)

cet_build_plugin(G4BeamlineGenerator art::module REG_SOURCE src/G4BeamlineGenerator_module.cc)

cet_build_plugin(GeneratorPlots art::module REG_SOURCE src/GeneratorPlots_module.cc)

cet_build_plugin(GenEventCounter art::module REG_SOURCE src/GenEventCounter_module.cc)

cet_build_plugin(InFlightParticleSampler art::module REG_SOURCE src/InFlightParticleSampler_module.cc)

cet_build_plugin(LeadingLog art::module REG_SOURCE src/LeadingLog_module.cc)

cet_build_plugin(MuStopProductsGun art::module REG_SOURCE src/MuStopProductsGun_module.cc)

cet_build_plugin(PBIWeight art::module REG_SOURCE src/PBIWeight_module.cc)

cet_build_plugin(Pileup art::module REG_SOURCE src/Pileup_module.cc)

cet_build_plugin(PrimaryAntiProtonGun art::module REG_SOURCE src/PrimaryAntiProtonGun_module.cc)

cet_build_plugin(PrimaryProtonGun art::module REG_SOURCE src/PrimaryProtonGun_module.cc)

cet_build_plugin(RanTest art::module REG_SOURCE src/RanTest_module.cc)

cet_build_plugin(RMCGun art::module REG_SOURCE src/RMCGun_module.cc)

cet_build_plugin(RPCGun art::module REG_SOURCE src/RPCGun_module.cc)

cet_build_plugin(SingleProcessGenerator art::module REG_SOURCE src/SingleProcessGenerator_module.cc)

cet_build_plugin(StoppedMuonRMCGun art::module REG_SOURCE src/StoppedMuonRMCGun_module.cc)

cet_build_plugin(StoppedMuonXRayGammaRayGun art::module REG_SOURCE src/StoppedMuonXRayGammaRayGun_module.cc)

cet_build_plugin(StoppedMuplusDecayGun art::module REG_SOURCE src/StoppedMuplusDecayGun_module.cc)

cet_build_plugin(StoppedParticleG4Gun art::module REG_SOURCE src/StoppedParticleG4Gun_module.cc)

cet_build_plugin(StoppedParticleReactionGun art::module REG_SOURCE src/StoppedParticleReactionGun_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME EventGenerator)
