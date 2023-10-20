cet_make_library(
    SOURCE
      src/CompressionLevel.cc
    LIBRARIES PUBLIC
      
)

cet_build_plugin(CompressDetStepMCs art::module
    REG_SOURCE src/CompressDetStepMCs_module.cc
    LIBRARIES REG
)

cet_build_plugin(CompressDigiMCsCheck art::module
    REG_SOURCE src/CompressDigiMCsCheck_module.cc
    LIBRARIES REG
)

cet_build_plugin(CompressDigiMCs art::module
    REG_SOURCE src/CompressDigiMCs_module.cc
    LIBRARIES REG
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Compression)
