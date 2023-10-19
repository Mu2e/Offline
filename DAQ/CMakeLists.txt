cet_make_library(
    SOURCE
src/CaloDAQUtilities.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(ArtBinaryPacketsFromDigis art::module REG_SOURCE src/ArtBinaryPacketsFromDigis_module.cc)

cet_build_plugin(ArtFragmentsFromDTCEvents art::module REG_SOURCE src/ArtFragmentsFromDTCEvents_module.cc)

cet_build_plugin(CaloHitsFromFragments art::module REG_SOURCE src/CaloHitsFromFragments_module.cc)

cet_build_plugin(CaloRecoFromFragments art::module REG_SOURCE src/CaloRecoFromFragments_module.cc)

cet_build_plugin(CrvDigisFromFragments art::module REG_SOURCE src/CrvDigisFromFragments_module.cc)

cet_build_plugin(FragmentAna art::module REG_SOURCE src/FragmentAna_module.cc)

cet_build_plugin(PrefetchDAQData art::module REG_SOURCE src/PrefetchDAQData_module.cc)

cet_build_plugin(StrawHitRecoFromFragments art::module REG_SOURCE src/StrawHitRecoFromFragments_module.cc)

cet_build_plugin(StrawRecoFromFragments art::module REG_SOURCE src/StrawRecoFromFragments_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME DAQ)
