cet_make_library(
    SOURCE
      src/CompressionLevel.cc
    LIBRARIES PUBLIC
      Boost::headers
      
)

cet_build_plugin(CompressDetStepMCs art::module
    REG_SOURCE src/CompressDetStepMCs_module.cc
    LIBRARIES REG
      Offline::Compression
      
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

cet_build_plugin(CompressDigiMCsCheck art::module
    REG_SOURCE src/CompressDigiMCsCheck_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::Compression
      
      Offline::DataProducts
      Offline::MCDataProducts
)

cet_build_plugin(CompressDigiMCs art::module
    REG_SOURCE src/CompressDigiMCs_module.cc
    LIBRARIES REG
      Offline::Compression
      
      Offline::DataProducts
      Offline::MCDataProducts
      Offline::Mu2eUtilities
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/Compression/fcl)
