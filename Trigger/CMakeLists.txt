cet_build_plugin(DigiFilter art::module
    REG_SOURCE src/DigiFilter_module.cc
    LIBRARIES REG
      Offline::RecoDataProducts
)

cet_build_plugin(EvalWeightedEventCounts art::module
    REG_SOURCE src/EvalWeightedEventCounts_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
)

cet_build_plugin(MergeTriggerInfo art::module
    REG_SOURCE src/MergeTriggerInfo_module.cc
    LIBRARIES REG
      Offline::RecoDataProducts
)

cet_build_plugin(PrescaleEvent art::module
    REG_SOURCE src/PrescaleEvent_module.cc
    LIBRARIES REG
      artdaq-core-mu2e::Data
)

cet_build_plugin(ReadTriggerInfo art::module
    REG_SOURCE src/ReadTriggerInfo_module.cc
    LIBRARIES REG
      
      Offline::BFieldGeom
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(ReadTriggerPath art::module
    REG_SOURCE src/ReadTriggerPath_module.cc
    LIBRARIES REG
      
      Offline::BFieldGeom
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

install_source(SUBDIRS src)
