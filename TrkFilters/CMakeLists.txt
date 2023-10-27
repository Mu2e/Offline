cet_build_plugin(HelixFilter art::module
    REG_SOURCE src/HelixFilter_module.cc
    LIBRARIES REG      
      Offline::BFieldGeom
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

cet_build_plugin(SeedFilter art::module
    REG_SOURCE src/SeedFilter_module.cc
    LIBRARIES REG      
      Offline::RecoDataProducts
)

cet_build_plugin(TimeClusterFilter art::module
    REG_SOURCE src/TimeClusterFilter_module.cc
    LIBRARIES REG      
      Offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrkFilters)
