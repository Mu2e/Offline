cet_make_library(
    SOURCE
src/TrkExtDetectors.cc
src/TrkExtDiag.cc
src/TrkExtInstanceName.cc
src/TrkExtMaterial.cc
src/TrkExtMCHits.cc
src/TrkExtProtonAbsorber.cc
src/TrkExtShape.cc
src/TrkExtStoppingTarget.cc
src/TrkExtToyDS.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(TrkExt art::module REG_SOURCE src/TrkExt_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
