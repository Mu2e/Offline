cet_build_plugin(HelixFilter art::module
    REG_SOURCE src/HelixFilter_module.cc
    LIBRARIES REG
)

cet_build_plugin(SeedFilter art::module
    REG_SOURCE src/SeedFilter_module.cc
    LIBRARIES REG
)

cet_build_plugin(TimeClusterFilter art::module
    REG_SOURCE src/TimeClusterFilter_module.cc
    LIBRARIES REG
)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME TrkFilters)
