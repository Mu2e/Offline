cet_make_library(
    SOURCE
      src/CrvHelper.cc
      src/MakeCrvRecoPulses.cc
    LIBRARIES PUBLIC
      
      offline::CosmicRayShieldGeom
      offline::DataProducts
      offline::GeometryService
      ROOT::Hist
)

cet_build_plugin(CrvCoincidenceFinder art::module
    REG_SOURCE src/CrvCoincidenceFinder_module.cc
    LIBRARIES REG
      offline::CRVReco
      
      offline::CosmicRayShieldGeom
      offline::DataProducts
      offline::GeometryService
      offline::RecoDataProducts
)

cet_build_plugin(CrvRecoPulsesFinder art::module
    REG_SOURCE src/CrvRecoPulsesFinder_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::CRVReco
      
      offline::CosmicRayShieldGeom
      offline::CRVConditions
      offline::DataProducts
      offline::GeometryService
      offline::ProditionsService
      offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CRVReco)
