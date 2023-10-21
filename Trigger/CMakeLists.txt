cet_build_plugin(DigiFilter art::module
    REG_SOURCE src/DigiFilter_module.cc
    LIBRARIES REG
      offline::RecoDataProducts
)

cet_build_plugin(EvalWeightedEventCounts art::module
    REG_SOURCE src/EvalWeightedEventCounts_module.cc
    LIBRARIES REG
      offline::MCDataProducts
)

cet_build_plugin(MergeTriggerInfo art::module
    REG_SOURCE src/MergeTriggerInfo_module.cc
    LIBRARIES REG
      offline::RecoDataProducts
)

cet_build_plugin(PrescaleEvent art::module
    REG_SOURCE src/PrescaleEvent_module.cc
    LIBRARIES REG
      
)

cet_build_plugin(ReadTriggerInfo art::module
    REG_SOURCE src/ReadTriggerInfo_module.cc
    LIBRARIES REG
      
      offline::BFieldGeom
      offline::ConditionsService
      offline::DataProducts
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::TrackerGeom
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Trigger)
