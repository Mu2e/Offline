cet_build_plugin(CountMixedEvents art::module
    REG_SOURCE src/CountMixedEvents_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(DetectorEfficiency art::module
    REG_SOURCE src/HPGeEfficiency_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(HPGeTree_module art::module
    REG_SOURCE src/HPGeTree_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(HPGeWaveformsFromStepPointMCs art::module
    REG_SOURCE src/HPGeWaveformsFromStepPointMCs_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::STMConditions
      Offline::GlobalConstantsService
)

cet_build_plugin(ShiftVirtualDetectorStepPointMCs art::module
    REG_SOURCE src/ShiftVirtualDetectorStepPointMCs_module.cc
    LIBRARIES REG
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

cet_build_plugin(STMResamplingFilter art::module
    REG_SOURCE src/STMResamplingProducer_module.cc
    LIBRARIES REG
        Offline::MCDataProducts
)

cet_build_plugin(STMResamplingProducer art::module
    REG_SOURCE src/STMResamplingProducer_module.cc
    LIBRARIES REG
        Offline::MCDataProducts
)

cet_build_plugin(VirtualDetectorTree_module art::module
    REG_SOURCE src/VirtualDetectorTree_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/HPGeWaveformsFromStepPointMCs.fcl   ${CURRENT_BINARY_DIR} fcl/HPGeWaveformsFromStepPointMCs.fcl   COPYONLY)


install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/STMMC/fcl)
