cet_build_plugin(ConcatenateDigitizedWaveforms art::module
    REG_SOURCE src/ConcatenateDigitizedWaveforms_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::RecoDataProducts
)

cet_build_plugin(ConcatenationFilter art::module
    REG_SOURCE src/ConcatenationFilter_module.cc
    LIBRARIES REG
    Offline::RecoDataProducts
)

cet_build_plugin(CountMixedEvents art::module
    REG_SOURCE src/CountMixedEvents_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(DetectorEfficiency art::module
    REG_SOURCE src/DetectorEfficiency_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(HPGeTree art::module
    REG_SOURCE src/HPGeTree_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(HPGeWaveformsFromStepPointMCs art::module
    REG_SOURCE src/HPGeWaveformsFromStepPointMCs_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::STMConditions
      Offline::GlobalConstantsService
)
cet_build_plugin(MWDTree art::module
    REG_SOURCE src/MWDTree_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::RecoDataProducts
)

cet_build_plugin(ShiftVirtualDetectorStepPointMCs art::module
    REG_SOURCE src/ShiftVirtualDetectorStepPointMCs_module.cc
    LIBRARIES REG
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)

cet_build_plugin(STMResamplingProducer art::module
    REG_SOURCE src/STMResamplingProducer_module.cc
    LIBRARIES REG
        Offline::MCDataProducts
)

cet_build_plugin(STMResamplingFilter art::module
    REG_SOURCE src/STMResamplingFilter_module.cc
    LIBRARIES REG
        Offline::MCDataProducts
)

cet_build_plugin(ValidateConsecutiveEventIDs art::module
    REG_SOURCE src/ValidateConsecutiveEventIDs_module.cc
    LIBRARIES REG
      Offline::MCDataProducts
)

cet_build_plugin(VirtualDetectorTree art::module
    REG_SOURCE src/VirtualDetectorTree_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GlobalConstantsService
      Offline::MCDataProducts
)


install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/STMMC/fcl)
