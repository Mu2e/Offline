cet_make_library(
    SOURCE
src/CaloRawWFProcessor.cc
src/CaloTemplateWFProcessor.cc
src/CaloTemplateWFUtil.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(CaloHitMakerFast art::module REG_SOURCE src/CaloHitMakerFast_module.cc)

cet_build_plugin(CaloHitMaker art::module REG_SOURCE src/CaloHitMaker_module.cc)

cet_build_plugin(CaloRecoDigiMaker art::module REG_SOURCE src/CaloRecoDigiMaker_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CaloReco)
