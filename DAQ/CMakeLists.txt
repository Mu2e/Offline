cet_make_library(
    SOURCE
      src/CaloDAQUtilities.cc
    LIBRARIES PUBLIC
      artdaq-core-mu2e::Data

)

cet_build_plugin(ArtBinaryPacketsFromDigis art::module
    REG_SOURCE src/ArtBinaryPacketsFromDigis_module.cc
    LIBRARIES REG
      Offline::DAQ

      Offline::CaloConditions
      Offline::CalorimeterGeom
      Offline::CosmicRayShieldGeom
      Offline::CRVConditions
      Offline::GeometryService
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::SeedService
)

cet_build_plugin(ProcessDTCAndCFOEvents art::module
    REG_SOURCE src/ProcessDTCAndCFOEvents_module.cc
    LIBRARIES REG
      Offline::DAQ
)

cet_build_plugin(EventHeaderFromCFOFragment art::module
    REG_SOURCE src/EventHeaderFromCFOFragment_module.cc
    LIBRARIES REG
      Offline::DAQ
)
    

cet_build_plugin(CaloHitsFromDataDecoders art::module
    REG_SOURCE src/CaloHitsFromDataDecoders_module.cc
    LIBRARIES REG
      Offline::DAQ

      Offline::CaloConditions
      Offline::ProditionsService
      Offline::RecoDataProducts
)

cet_build_plugin(CaloDigisFromDataDecoders art::module
    REG_SOURCE src/CaloDigisFromDataDecoders_module.cc
    LIBRARIES REG
      Offline::DAQ

      Offline::CaloConditions
      Offline::ProditionsService
      Offline::RecoDataProducts
      artdaq-core-mu2e::Data
      artdaq-core-mu2e::Data_dict
)

cet_build_plugin(CRVDigisFromDataDecoders art::module
    REG_SOURCE src/CRVDigisFromDataDecoders_module.cc
    LIBRARIES REG
      Offline::DAQ

      Offline::CRVConditions
      Offline::ProditionsService
      Offline::RecoDataProducts
      artdaq-core-mu2e::Data
      artdaq-core-mu2e::Data_dict
)

cet_build_plugin(CRVGlobalRunDataFromDataDecoders art::module
    REG_SOURCE src/CRVGlobalRunDataFromDataDecoders_module.cc
    LIBRARIES REG
      Offline::DAQ

      Offline::CRVConditions
      Offline::ProditionsService
      Offline::RecoDataProducts
      artdaq-core-mu2e::Data
      artdaq-core-mu2e::Data_dict
)

cet_build_plugin(DataDecodersAna art::module
    REG_SOURCE src/DataDecodersAna_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DAQ
      Offline::DataProducts
      ROOT::Hist
)

cet_build_plugin(DummyLumiInfoProducer art::module
    REG_SOURCE src/DummyLumiInfoProducer_module.cc
    LIBRARIES REG
      Offline::DAQ
      Offline::RecoDataProducts
)

cet_build_plugin(LumiInfoAna art::module
    REG_SOURCE src/LumiInfoAna_module.cc
    LIBRARIES REG
    Offline::DAQ
    Offline::RecoDataProducts
    art_root_io::TFileService_service
    artdaq-core-mu2e::Data
    artdaq-core-mu2e::Data_dict
)

cet_build_plugin(PrefetchDAQData art::module
    REG_SOURCE src/PrefetchDAQData_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DAQ

      Offline::RecoDataProducts
)

cet_build_plugin(ComboHitsFromDataDecoders art::module
    REG_SOURCE src/ComboHitsFromDataDecoders_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DAQ

      Offline::ConditionsBase
      Offline::ConditionsService
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeometryService
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::TrackerConditions
      Offline::TrackerGeom
      Offline::TrkHitReco
)

cet_build_plugin(StrawDigisFromDataDecoders art::module
    REG_SOURCE src/StrawDigisFromDataDecoders_module.cc
    LIBRARIES REG
      Offline::DAQ

      Offline::DataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(STMWaveformDigisFromFragments art::module
    REG_SOURCE src/STMWaveformDigisFromFragments_module.cc
    LIBRARIES REG
      Offline::DAQ
      Offline::ProditionsService
      Offline::DataProducts
      Offline::RecoDataProducts
      BTrk_difAlgebra
)

cet_build_plugin(StrawDigisFromArtdaqFragments art::module
  REG_SOURCE src/StrawDigisFromArtdaqFragments_module.cc
  LIBRARIES REG
  Offline::DAQ
  Offline::DataProducts
  Offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/DAQ/fcl)
