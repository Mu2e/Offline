cet_build_plugin(CrvCoincidenceClusterFilter art::module
    REG_SOURCE src/CrvCoincidenceClusterFilter_module.cc
    LIBRARIES REG
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME CRVFilters)
