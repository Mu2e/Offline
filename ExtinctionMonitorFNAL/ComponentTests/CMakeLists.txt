#cet_build_plugin(PixelIdConverterTest art::module
#    REG_SOURCE src/PixelIdConverterTest_module.cc
#    LIBRARIES REG
#      Offline::DataProducts
#      Offline::GeometryService
#      Offline::SeedService
#)

#cet_build_plugin(TrackExtrapolatorTest2 art::module
#    REG_SOURCE src/TrackExtrapolatorTest2_module.cc
#    LIBRARIES REG
#      Offline::GeometryService
#      Offline::MCDataProducts
#      Offline::RecoDataProducts
#      art_root_io::TFileService_service
#)

#cet_build_plugin(TrackExtrapolatorTest art::module
#    REG_SOURCE src/TrackExtrapolatorTest_module.cc
#    LIBRARIES REG
#      Offline::GeometryService
#      Offline::RecoDataProducts
#)

install_source(SUBDIRS src)
