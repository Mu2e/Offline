cet_build_plugin(PixelIdConverterTest art::module
    REG_SOURCE src/PixelIdConverterTest_module.cc
    LIBRARIES REG
      offline::DataProducts
      offline::GeometryService
      offline::SeedService
)

cet_build_plugin(TrackExtrapolatorTest2 art::module
    REG_SOURCE src/TrackExtrapolatorTest2_module.cc
    LIBRARIES REG
      offline::GeometryService
      offline::MCDataProducts
      offline::RecoDataProducts
      art_root_io::TFileService_service
)

cet_build_plugin(TrackExtrapolatorTest art::module
    REG_SOURCE src/TrackExtrapolatorTest_module.cc
    LIBRARIES REG
      offline::GeometryService
      offline::RecoDataProducts
)

install_source(SUBDIRS src)
