cet_make_library(
    SOURCE
      src/CaloDAQUtilities.cc
    LIBRARIES PUBLIC
      artdaq_core_mu2e::artdaq-core-mu2e_Data
      
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

cet_build_plugin(ArtFragmentsFromDTCEvents art::module
    REG_SOURCE src/ArtFragmentsFromDTCEvents_module.cc
    LIBRARIES REG
      Offline::DAQ
      
)

cet_build_plugin(CaloHitsFromFragments art::module
    REG_SOURCE src/CaloHitsFromFragments_module.cc
    LIBRARIES REG
      Offline::DAQ
      
      Offline::CaloConditions
      Offline::ProditionsService
      Offline::RecoDataProducts
)

cet_build_plugin(CaloDigiFromFragments art::module
    REG_SOURCE src/CaloDigiFromFragments_module.cc
    LIBRARIES REG
      Offline::DAQ
      
      Offline::CaloConditions
      Offline::ProditionsService
      Offline::RecoDataProducts
)

cet_build_plugin(CrvDigisFromFragments art::module
    REG_SOURCE src/CrvDigisFromFragments_module.cc
    LIBRARIES REG
      Offline::DAQ
      
      Offline::CRVConditions
      Offline::ProditionsService
      Offline::RecoDataProducts
)

cet_build_plugin(FragmentAna art::module
    REG_SOURCE src/FragmentAna_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DAQ
      Offline::DataProducts
      ROOT::Hist
)

cet_build_plugin(PrefetchDAQData art::module
    REG_SOURCE src/PrefetchDAQData_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::DAQ
      
      Offline::RecoDataProducts
)

cet_build_plugin(StrawHitRecoFromFragments art::module
    REG_SOURCE src/StrawHitRecoFromFragments_module.cc
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

cet_build_plugin(StrawRecoFromFragments art::module
    REG_SOURCE src/StrawRecoFromFragments_module.cc
    LIBRARIES REG
      Offline::DAQ
      
      Offline::DataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(StrawRecoFromArtdaqFragments art::module
  REG_SOURCE src/StrawRecoFromArtdaqFragments_module.cc
  LIBRARIES REG
  Offline::DAQ
  Offline::DataProducts
  Offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/DAQ/fcl)
