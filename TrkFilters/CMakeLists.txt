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

cet_build_plugin(KalSeedFilter art::module
    REG_SOURCE src/KalSeedFilter_module.cc
    LIBRARIES REG      
      Offline::RecoDataProducts
)

cet_build_plugin(MultiHelixFilter art::module
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

cet_build_plugin(TimeClusterFilter art::module
    REG_SOURCE src/TimeClusterFilter_module.cc
    LIBRARIES REG      
      Offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/TrkFilters/fcl)
