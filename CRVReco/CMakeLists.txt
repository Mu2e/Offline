cet_make_library(
    SOURCE
      src/CrvHelper.cc
      src/MakeCrvRecoPulses.cc
    LIBRARIES PUBLIC
      
      Offline::CosmicRayShieldGeom
      Offline::DataProducts
      Offline::GeometryService
      ROOT::Hist
)

cet_build_plugin(CrvCoincidenceFinder art::module
    REG_SOURCE src/CrvCoincidenceFinder_module.cc
    LIBRARIES REG
      Offline::CRVReco
      
      Offline::CosmicRayShieldGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(CrvRecoPulsesFinder art::module
    REG_SOURCE src/CrvRecoPulsesFinder_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CRVReco
      
      Offline::CosmicRayShieldGeom
      Offline::CRVConditions
      Offline::DataProducts
      Offline::GeometryService
      Offline::ProditionsService
      Offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CRVReco)
