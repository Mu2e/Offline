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
      
      art_root_io::tfile_support
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::DataProducts
      Offline::ExtinctionMonitorFNAL_Geometry
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProductionTargetGeom
      Offline::SeedService
      Offline::StoppingTargetGeom
      CRY
)

cet_build_plugin(AntiProtonGun art::module
    REG_SOURCE src/AntiProtonGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProductionTargetGeom
      Offline::SeedService
)

cet_build_plugin(AntiprotonResampling art::module
    REG_SOURCE src/AntiprotonResampling_module.cc
    LIBRARIES REG
      Offline::EventGenerator
)

cet_build_plugin(CaloCalibGun art::module
    REG_SOURCE src/CaloCalibGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::CalorimeterGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(CaloTBGun art::module
    REG_SOURCE src/CaloTBGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(CeEndpoint art::module
    REG_SOURCE src/CeEndpoint_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(CORSIKAEventGenerator art::module
    REG_SOURCE src/CORSIKAEventGenerator_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::ExtinctionMonitorFNAL_Geometry
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(CRYEventGenerator art::module
    REG_SOURCE src/CRYEventGenerator_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(CryResampler art::module
    REG_SOURCE src/CryResampler_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(DecayInOrbitWeight art::module
    REG_SOURCE src/DecayInOrbitWeight_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(deltaFunctionGun art::module
    REG_SOURCE src/deltaFunctionGun_module.cc
    LIBRARIES REG
        Offline::EventGenerator
)

cet_build_plugin(EplusFromStoppedPion art::module
    REG_SOURCE src/EplusFromStoppedPion_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(EventGenerator art::module
    REG_SOURCE src/EventGenerator_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(ExtMonFNALBoxGenerator art::module
    REG_SOURCE src/ExtMonFNALBoxGenerator_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::DataProducts
      Offline::ExtinctionMonitorFNAL_Geometry
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(ExtMonFNALGun art::module
    REG_SOURCE src/ExtMonFNALGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(ExtMonFNALRoomGenerator art::module
    REG_SOURCE src/ExtMonFNALRoomGenerator_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::DataProducts
      Offline::ExtinctionMonitorFNAL_Geometry
      Offline::ExtinctionMonitorFNAL_Utilities
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::ProtonBeamDumpGeom
      Offline::SeedService
)

cet_build_plugin(FlatGeneratorFromStoppedMuon art::module
    REG_SOURCE src/FlatGeneratorFromStoppedMuon_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(FromAsciiMomentumAndPosition art::module
    REG_SOURCE src/FromAsciiMomentumAndPosition_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProductionTargetGeom
      Offline::SeedService
)

cet_build_plugin(FromSimParticleCompact art::module
    REG_SOURCE src/FromSimParticleCompact_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProductionTargetGeom
      Offline::SeedService
)

cet_build_plugin(FromStepPointMCs art::module
    REG_SOURCE src/FromStepPointMCs_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

cet_build_plugin(FromStepPointMCsRotateTarget art::module
    REG_SOURCE src/FromStepPointMCsRotateTarget_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::ProductionTargetGeom
      Offline::SeedService
)

cet_build_plugin(G4BeamlineGenerator art::module
    REG_SOURCE src/G4BeamlineGenerator_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(GammaConvFlat art::module
    REG_SOURCE src/GammaConvFlat_module.cc
    LIBRARIES REG
      Offline::EventGenerator
)

cet_build_plugin(GammaConversion art::module
    REG_SOURCE src/GammaConversion_module.cc
    LIBRARIES REG
      Offline::EventGenerator
)

cet_build_plugin(GeneratorPlots art::module
    REG_SOURCE src/GeneratorPlots_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::MCDataProducts
)

cet_build_plugin(GenEventCounter art::module
    REG_SOURCE src/GenEventCounter_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::MCDataProducts
)

cet_build_plugin(GenFilter art::module
    REG_SOURCE src/GenFilter_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::MCDataProducts
)

cet_build_plugin(InFlightParticleSampler art::module
    REG_SOURCE src/InFlightParticleSampler_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(LeadingLog art::module
    REG_SOURCE src/LeadingLog_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(MuStopProductsGun art::module
    REG_SOURCE src/MuStopProductsGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(PBIWeight art::module
    REG_SOURCE src/PBIWeight_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::MCDataProducts
)

cet_build_plugin(Pileup art::module
    REG_SOURCE src/Pileup_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(PionFilter art::module
    REG_SOURCE src/PionFilter_module.cc
    LIBRARIES REG
      Offline::EventGenerator
)

cet_build_plugin(PrimaryAntiProtonGun art::module
    REG_SOURCE src/PrimaryAntiProtonGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProductionTargetGeom
      Offline::SeedService
)

cet_build_plugin(PrimaryProtonGun art::module
    REG_SOURCE src/PrimaryProtonGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(RanTest art::module
    REG_SOURCE src/RanTest_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::MCDataProducts
      Offline::SeedService
)

cet_build_plugin(RMCGun art::module
    REG_SOURCE src/RMCGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(RPCGun art::module
    REG_SOURCE src/RPCGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(SimpleAntiprotonGun art::module
    REG_SOURCE src/SimpleAntiprotonGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
)

cet_build_plugin(SingleProcessGenerator art::module
    REG_SOURCE src/SingleProcessGenerator_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(StoppedMuonRMCGun art::module
    REG_SOURCE src/StoppedMuonRMCGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(StoppedMuonXRayGammaRayGun art::module
    REG_SOURCE src/StoppedMuonXRayGammaRayGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
      Offline::StoppingTargetGeom
)

cet_build_plugin(StoppedMuplusDecayGun art::module
    REG_SOURCE src/StoppedMuplusDecayGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(StoppedParticleG4Gun art::module
    REG_SOURCE src/StoppedParticleG4Gun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(StoppedParticleReactionGun art::module
    REG_SOURCE src/StoppedParticleReactionGun_module.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::SeedService
)

cet_build_plugin(DIOGenerator art::tool
    REG_SOURCE src/DIOGenerator_tool.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(Mu2eXGenerator art::tool
    REG_SOURCE src/Mu2eXGenerator_tool.cc
    LIBRARIES REG
      Offline::EventGenerator

      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(MuCap1809keVGammaGenerator art::tool
    REG_SOURCE src/MuCap1809keVGammaGenerator_tool.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(MuCapDeuteronGenerator art::tool
    REG_SOURCE src/MuCapDeuteronGenerator_tool.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(MuCapNeutronGenerator art::tool
    REG_SOURCE src/MuCapNeutronGenerator_tool.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(MuCapPhotonGenerator art::tool
    REG_SOURCE src/MuCapPhotonGenerator_tool.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(MuCapProtonGenerator art::tool
    REG_SOURCE src/MuCapProtonGenerator_tool.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(MuStopDecayPositionSamplerTool art::tool
    REG_SOURCE src/MuStopDecayPositionSamplerTool_tool.cc
    LIBRARIES REG
      Offline::EventGenerator      

      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(MuplusMichelGenerator art::tool
    REG_SOURCE src/MuplusMichelGenerator_tool.cc
    LIBRARIES REG
      Offline::EventGenerator
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(PhotonGun art::module
    REG_SOURCE src/PhotonGun_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::DataProducts
)

cet_build_plugin(RMCGenerator art::tool
    REG_SOURCE src/RMCGenerator_tool.cc
    LIBRARIES REG
      Offline::EventGenerator
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defaultConfigs/cosmicDYB.txt	      ${CURRENT_BINARY_DIR} defaultConfigs/cosmicDYB.txt	     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defaultConfigs/defaultCRYconfig.txt   ${CURRENT_BINARY_DIR} defaultConfigs/defaultCRYconfig.txt  )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defaultConfigs/extractedCRYconfig.txt ${CURRENT_BINARY_DIR} defaultConfigs/extractedCRYconfig.txt)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defaultConfigs/fromG4BLFile.txt	      ${CURRENT_BINARY_DIR} defaultConfigs/fromG4BLFile.txt	     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defaultConfigs/muCRYconfig.txt	      ${CURRENT_BINARY_DIR} defaultConfigs/muCRYconfig.txt	     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defaultConfigs/othersCRYconfig.txt    ${CURRENT_BINARY_DIR} defaultConfigs/othersCRYconfig.txt   )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defaultConfigs/particleGun.txt	      ${CURRENT_BINARY_DIR} defaultConfigs/particleGun.txt	     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defaultConfigs/primaryProtonGun.txt   ${CURRENT_BINARY_DIR} defaultConfigs/primaryProtonGun.txt  )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/defaultConfigs/surfaceCheck.txt       ${CURRENT_BINARY_DIR} defaultConfigs/surfaceCheck.txt      )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/acDipoleTransmissionFunction.txt      ${CURRENT_BINARY_DIR} data/acDipoleTransmissionFunction.txt     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/czarnecki_Al.tbl		      ${CURRENT_BINARY_DIR} data/czarnecki_Al.tbl		     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/czarnecki_szafron_Al_2016_tail.tbl    ${CURRENT_BINARY_DIR} data/czarnecki_szafron_Al_2016_tail.tbl   )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/czarnecki_szafron_Al_2016.tbl	      ${CURRENT_BINARY_DIR} data/czarnecki_szafron_Al_2016.tbl	     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/czarnecki_Ti.tbl		      ${CURRENT_BINARY_DIR} data/czarnecki_Ti.tbl		     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/heeck_finer_binning_2016_szafron.tbl  ${CURRENT_BINARY_DIR} data/heeck_finer_binning_2016_szafron.tbl )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/Micheltbl.tbl			      ${CURRENT_BINARY_DIR} data/Micheltbl.tbl			     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/neutronSpectrum.txt		      ${CURRENT_BINARY_DIR} data/neutronSpectrum.txt		     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/neutSpectrumDocdb1619.txt	      ${CURRENT_BINARY_DIR} data/neutSpectrumDocdb1619.txt	     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/photonSpectrumMuonicAlCapture_v1.txt  ${CURRENT_BINARY_DIR} data/photonSpectrumMuonicAlCapture_v1.txt )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/PlettNeutronSpectrum.txt	      ${CURRENT_BINARY_DIR} data/PlettNeutronSpectrum.txt	     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/SchroederNeutronSpectrum.txt	      ${CURRENT_BINARY_DIR} data/SchroederNeutronSpectrum.txt	     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/StoppedMuons.txt		      ${CURRENT_BINARY_DIR} data/StoppedMuons.txt		     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/szafron_mu2e_7615.tbl		      ${CURRENT_BINARY_DIR} data/szafron_mu2e_7615.tbl		     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/timeDelayDist.txt                     ${CURRENT_BINARY_DIR} data/timeDelayDist.txt                   )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/watanabe.tbl                          ${CURRENT_BINARY_DIR} data/watanabe.tbl                        )

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/EventGenerator)
install(DIRECTORY defaultConfigs DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/EventGenerator)
install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/EventGenerator/fcl)
