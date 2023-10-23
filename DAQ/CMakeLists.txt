cet_make_library(
    SOURCE
      src/CaloDAQUtilities.cc
    LIBRARIES PUBLIC
      artdaq_core_mu2e::artdaq-core-mu2e_Data
      
)

cet_build_plugin(ArtBinaryPacketsFromDigis art::module
    REG_SOURCE src/ArtBinaryPacketsFromDigis_module.cc
    LIBRARIES REG
      offline::DAQ
      
      offline::CaloConditions
      offline::CalorimeterGeom
      offline::CosmicRayShieldGeom
      offline::CRVConditions
      offline::GeometryService
      offline::ProditionsService
      offline::RecoDataProducts
      offline::SeedService
)

cet_build_plugin(ArtFragmentsFromDTCEvents art::module
    REG_SOURCE src/ArtFragmentsFromDTCEvents_module.cc
    LIBRARIES REG
      offline::DAQ
      
)

cet_build_plugin(CaloHitsFromFragments art::module
    REG_SOURCE src/CaloHitsFromFragments_module.cc
    LIBRARIES REG
      offline::DAQ
      
      offline::CaloConditions
      offline::ProditionsService
      offline::RecoDataProducts
)

cet_build_plugin(CaloRecoFromFragments art::module
    REG_SOURCE src/CaloRecoFromFragments_module.cc
    LIBRARIES REG
      offline::DAQ
      
      offline::CaloConditions
      offline::ProditionsService
      offline::RecoDataProducts
)

cet_build_plugin(CrvDigisFromFragments art::module
    REG_SOURCE src/CrvDigisFromFragments_module.cc
    LIBRARIES REG
      offline::DAQ
      
      offline::CRVConditions
      offline::ProditionsService
      offline::RecoDataProducts
)

cet_build_plugin(FragmentAna art::module
    REG_SOURCE src/FragmentAna_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::DAQ
      offline::DataProducts
      ROOT::Hist
)

cet_build_plugin(PrefetchDAQData art::module
    REG_SOURCE src/PrefetchDAQData_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::DAQ
      
      offline::RecoDataProducts
)

cet_build_plugin(StrawHitRecoFromFragments art::module
    REG_SOURCE src/StrawHitRecoFromFragments_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::DAQ
      
      offline::ConditionsBase
      offline::ConditionsService
      offline::ConfigTools
      offline::DataProducts
      offline::GeometryService
      offline::ProditionsService
      offline::RecoDataProducts
      offline::TrackerConditions
      offline::TrackerGeom
      offline::TrkHitReco
)

cet_build_plugin(StrawRecoFromFragments art::module
    REG_SOURCE src/StrawRecoFromFragments_module.cc
    LIBRARIES REG
      offline::DAQ
      
      offline::DataProducts
      offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME DAQ)
