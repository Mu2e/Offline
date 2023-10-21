cet_make_library(
    SOURCE
      src/CosmicCRY.cc
      src/CosmicDYB.cc
      src/CosmicFromTH2.cc
      src/DYBGenerator.cc
      src/ExtMonFNALGunImpl.cc
      src/FromG4BLFile.cc
      src/ParticleGun.cc
      src/ParticleGunImpl.cc
      src/PiCaptureEffects.cc
      src/PrimaryProtonGunImpl.cc
    LIBRARIES PUBLIC
      
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::ConfigTools
      offline::DataProducts
      offline::ExtinctionMonitorFNAL_Geometry
      offline::GeneralUtilities
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::ProductionTargetGeom
      offline::SeedService
      offline::StoppingTargetGeom
)

cet_build_plugin(AntiProtonGun art::module
    REG_SOURCE src/AntiProtonGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::DataProducts
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::ProductionTargetGeom
      offline::SeedService
)

cet_build_plugin(CaloCalibGun art::module
    REG_SOURCE src/CaloCalibGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::CalorimeterGeom
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(CaloTBGun art::module
    REG_SOURCE src/CaloTBGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::CalorimeterGeom
      offline::ConfigTools
      offline::DataProducts
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::SeedService
)

cet_build_plugin(CeEndpoint art::module
    REG_SOURCE src/CeEndpoint_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(CORSIKAEventGenerator art::module
    REG_SOURCE src/CORSIKAEventGenerator_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::CalorimeterGeom
      offline::ConfigTools
      offline::ExtinctionMonitorFNAL_Geometry
      offline::GeneralUtilities
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(CRYEventGenerator art::module
    REG_SOURCE src/CRYEventGenerator_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::MCDataProducts
      offline::SeedService
)

cet_build_plugin(CryResampler art::module
    REG_SOURCE src/CryResampler_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::GeneralUtilities
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::SeedService
)

cet_build_plugin(DecayInOrbitWeight art::module
    REG_SOURCE src/DecayInOrbitWeight_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(EplusFromStoppedPion art::module
    REG_SOURCE src/EplusFromStoppedPion_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::GeneralUtilities
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(EventGenerator art::module
    REG_SOURCE src/EventGenerator_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::MCDataProducts
      offline::SeedService
)

cet_build_plugin(ExtMonFNALBoxGenerator art::module
    REG_SOURCE src/ExtMonFNALBoxGenerator_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConditionsService
      offline::ConfigTools
      offline::DataProducts
      offline::ExtinctionMonitorFNAL_Geometry
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::SeedService
)

cet_build_plugin(ExtMonFNALGun art::module
    REG_SOURCE src/ExtMonFNALGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::MCDataProducts
      offline::SeedService
)

cet_build_plugin(ExtMonFNALRoomGenerator art::module
    REG_SOURCE src/ExtMonFNALRoomGenerator_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::DataProducts
      offline::ExtinctionMonitorFNAL_Geometry
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::ProtonBeamDumpGeom
      offline::SeedService
)

cet_build_plugin(FlatMuonDaughterGenerator art::module
    REG_SOURCE src/FlatMuonDaughterGenerator_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(FromAsciiMomentumAndPosition art::module
    REG_SOURCE src/FromAsciiMomentumAndPosition_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConditionsService
      offline::ConfigTools
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::ProductionTargetGeom
      offline::SeedService
)

cet_build_plugin(FromSimParticleCompact art::module
    REG_SOURCE src/FromSimParticleCompact_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::DataProducts
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::ProductionTargetGeom
      offline::SeedService
)

cet_build_plugin(FromStepPointMCs art::module
    REG_SOURCE src/FromStepPointMCs_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::GlobalConstantsService
      offline::MCDataProducts
)

cet_build_plugin(FromStepPointMCsRotateTarget art::module
    REG_SOURCE src/FromStepPointMCsRotateTarget_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::ProductionTargetGeom
      offline::SeedService
)

cet_build_plugin(G4BeamlineGenerator art::module
    REG_SOURCE src/G4BeamlineGenerator_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::MCDataProducts
      offline::SeedService
)

cet_build_plugin(GeneratorPlots art::module
    REG_SOURCE src/GeneratorPlots_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GeneralUtilities
      offline::MCDataProducts
)

cet_build_plugin(GenEventCounter art::module
    REG_SOURCE src/GenEventCounter_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::MCDataProducts
)

cet_build_plugin(InFlightParticleSampler art::module
    REG_SOURCE src/InFlightParticleSampler_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::GeneralUtilities
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(LeadingLog art::module
    REG_SOURCE src/LeadingLog_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(MuStopProductsGun art::module
    REG_SOURCE src/MuStopProductsGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::DataProducts
      offline::GeneralUtilities
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(PBIWeight art::module
    REG_SOURCE src/PBIWeight_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::MCDataProducts
)

cet_build_plugin(Pileup art::module
    REG_SOURCE src/Pileup_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(PrimaryAntiProtonGun art::module
    REG_SOURCE src/PrimaryAntiProtonGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConditionsService
      offline::ConfigTools
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::ProductionTargetGeom
      offline::SeedService
)

cet_build_plugin(PrimaryProtonGun art::module
    REG_SOURCE src/PrimaryProtonGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::MCDataProducts
      offline::SeedService
)

cet_build_plugin(RanTest art::module
    REG_SOURCE src/RanTest_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::MCDataProducts
      offline::SeedService
)

cet_build_plugin(RMCGun art::module
    REG_SOURCE src/RMCGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::DataProducts
      offline::GeneralUtilities
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(RPCGun art::module
    REG_SOURCE src/RPCGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(SingleProcessGenerator art::module
    REG_SOURCE src/SingleProcessGenerator_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(StoppedMuonRMCGun art::module
    REG_SOURCE src/StoppedMuonRMCGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::DataProducts
      offline::GeneralUtilities
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(StoppedMuonXRayGammaRayGun art::module
    REG_SOURCE src/StoppedMuonXRayGammaRayGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::DataProducts
      offline::GeneralUtilities
      offline::GeometryService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
      offline::StoppingTargetGeom
)

cet_build_plugin(StoppedMuplusDecayGun art::module
    REG_SOURCE src/StoppedMuplusDecayGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GeneralUtilities
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(StoppedParticleG4Gun art::module
    REG_SOURCE src/StoppedParticleG4Gun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GeneralUtilities
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(StoppedParticleReactionGun art::module
    REG_SOURCE src/StoppedParticleReactionGun_module.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::ConfigTools
      offline::DataProducts
      offline::GeneralUtilities
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::SeedService
)

cet_build_plugin(DIOGenerator art::tool
    REG_SOURCE src/DIOGenerator_tool.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(MuCap1809keVGammaGenerator art::tool
    REG_SOURCE src/MuCap1809keVGammaGenerator_tool.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(MuCapDeuteronGenerator art::tool
    REG_SOURCE src/MuCapDeuteronGenerator_tool.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(MuCapNeutronGenerator art::tool
    REG_SOURCE src/MuCapNeutronGenerator_tool.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(MuCapPhotonGenerator art::tool
    REG_SOURCE src/MuCapPhotonGenerator_tool.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(MuCapProtonGenerator art::tool
    REG_SOURCE src/MuCapProtonGenerator_tool.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
)

cet_build_plugin(MuplusMichelGenerator art::tool
    REG_SOURCE src/MuplusMichelGenerator_tool.cc
    LIBRARIES REG
      offline::EventGenerator
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME EventGenerator)
