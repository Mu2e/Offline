cet_build_plugin(TrackCaloIntersection art::module
    REG_SOURCE src/TrackCaloIntersection_module.cc
    LIBRARIES REG      
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(TrackCaloIntersectionMVA art::module
    REG_SOURCE src/TrackCaloIntersectionMVA_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(TrackCaloMatchingExtend art::module
    REG_SOURCE src/TrackCaloMatchingExtend_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(TrackCaloMatching art::module
    REG_SOURCE src/TrackCaloMatching_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::RecoDataProducts
      ROOT::Physics
)

cet_build_plugin(TrackCaloMatchingMVA art::module
    REG_SOURCE src/TrackCaloMatchingMVA_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(TrkExtrapol art::module
    REG_SOURCE src/TrkExtrapol_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
      Offline::TrackerGeom
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/TrackCaloMatching/fcl)
