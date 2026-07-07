cet_build_plugin(CrvDigiDump art::module
    REG_SOURCE src/CrvDigiDump_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers()
install_fhicl()

