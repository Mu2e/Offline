cet_make_library(
    SOURCE
      src/ClusterOnTrackPrecisionTool.cc
      src/LinearRegression.cc
      src/PixelHitLookup.cc
      src/PixelRecoUtils.cc
      src/TrackExtrapolator.cc
    LIBRARIES PUBLIC
      Offline::RecoDataProducts
)

cet_build_plugin(EMFPatRecFromTracklets art::module
    REG_SOURCE src/EMFPatRecFromTracklets_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ExtinctionMonitorFNAL_Reconstruction
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(EMFTrackClusterArbiter art::module
    REG_SOURCE src/EMFTrackClusterArbiter_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ExtinctionMonitorFNAL_Reconstruction
      
      Offline::RecoDataProducts
)

cet_build_plugin(ExtMonFNALRawClusterization art::module
    REG_SOURCE src/ExtMonFNALRawClusterization_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ExtinctionMonitorFNAL_Reconstruction
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(ExtMonFNALRecoClusterization art::module
    REG_SOURCE src/ExtMonFNALRecoClusterization_module.cc
    LIBRARIES REG
      Offline::ExtinctionMonitorFNAL_Reconstruction
      Offline::GeometryService
      Offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
