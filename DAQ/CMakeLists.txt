
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

cet_build_plugin(EventHeaderFromCFOFragment art::module
    REG_SOURCE src/EventHeaderFromCFOFragment_module.cc
    LIBRARIES REG
      Offline::DAQ
)


cet_build_plugin(CaloHitsFromDTCEvents art::module
    REG_SOURCE src/CaloHitsFromDTCEvents_module.cc
    LIBRARIES REG
      Offline::DAQ

      Offline::CaloConditions
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::DataProducts
)

cet_build_plugin(CaloDigisFromDTCEvents art::module
    REG_SOURCE src/CaloDigisFromDTCEvents_module.cc
    LIBRARIES REG
      Offline::DAQ

      Offline::CaloConditions
      Offline::ProditionsService
      Offline::RecoDataProducts
      artdaq-core-mu2e::Data
      artdaq-core-mu2e::Data_dict
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

cet_build_plugin(LumiStreamFilter art::module
  REG_SOURCE src/LumiStreamFilter_module.cc
  LIBRARIES REG
  fhiclcpp::types
  Offline::DAQ
  Offline::RecoDataProducts
  artdaq-core-mu2e::Data
  artdaq-core-mu2e::Data_dict
  TRACE::MF
  ${ARTDAQ_DAQDATA_LIB}
)

cet_build_plugin(PrefetchDAQData art::module
    REG_SOURCE src/PrefetchDAQData_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DAQ

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

cet_build_plugin(CrvDigisFromArtdaqFragments art::module
  REG_SOURCE src/CrvDigisFromArtdaqFragments_module.cc
  LIBRARIES REG
  Offline::DAQ
  Offline::DataProducts
  Offline::RecoDataProducts
  artdaq-core-mu2e::Data
  artdaq-core-mu2e::Overlays
)

cet_build_plugin(CrvGRdataFromArtdaqFragments art::module
  REG_SOURCE src/CrvGRdataFromArtdaqFragments_module.cc
  LIBRARIES REG
  Offline::DAQ
  Offline::DataProducts
  Offline::RecoDataProducts
  artdaq-core-mu2e::Data
  artdaq-core-mu2e::Overlays
)

cet_build_plugin(MTPHitsFromDTCEvents art::module
    REG_SOURCE src/MTPHitsFromDTCEvents_module.cc
    LIBRARIES REG
      Offline::DAQ
      Offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/DAQ/fcl)
