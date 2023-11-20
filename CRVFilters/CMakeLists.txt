cet_build_plugin(CrvCoincidenceClusterFilter art::module
    REG_SOURCE src/CrvCoincidenceClusterFilter_module.cc
    LIBRARIES REG
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CRVFilters/fcl)
