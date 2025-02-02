cet_build_plugin(HPGeWaveformsFromGeantSim art::module
    REG_SOURCE src/HPGeWaveformsFromGeantSim_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::DataProducts
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::STMConditions
      Offline::GlobalConstantsService
)
cet_build_plugin(STMResamplingProducer art::module
    REG_SOURCE src/STMResamplingProducer_module.cc
    LIBRARIES REG
        Offline::MCDataProducts
)
cet_build_plugin(STMResamplingFilter art::module
    REG_SOURCE src/STMResamplingProducer_module.cc
    LIBRARIES REG
        Offline::MCDataProducts
)
cet_build_plugin(DetectorEfficiency art::module
    REG_SOURCE src/HPGeEfficiency_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)
cet_build_plugin(EDepFromShower_module art::module
    REG_SOURCE src/EDepFromShower_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)
cet_build_plugin(VirtualDetectorTree_module art::module
    REG_SOURCE src/VirtualDetectorTree_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)
cet_build_plugin(ShiftVirtualDetectorStepPointMCs art::module
    REG_SOURCE src/ShiftVirtualDetectorStepPointMCs_module.cc
    LIBRARIES REG
      Offline::ConditionsService
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/HPGeWaveformsFromGeantSim.fcl   ${CURRENT_BINARY_DIR} fcl/HPGeWaveformsFromGeantSim.fcl   COPYONLY)


install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/STMMC/fcl)
