cet_make_library(
    SOURCE
      src/ClusterAssociator.cc
      src/ClusterFinder.cc
      src/ClusterUtils.cc
    LIBRARIES PUBLIC
      
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(CaloClusterFast art::module
    REG_SOURCE src/CaloClusterFast_module.cc
    LIBRARIES REG
      Offline::CaloCluster
      
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(CaloClusterMaker art::module
    REG_SOURCE src/CaloClusterMaker_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CaloCluster
      
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(CaloProtoClusterMaker art::module
    REG_SOURCE src/CaloProtoClusterMaker_module.cc
    LIBRARIES REG
      Offline::CaloCluster
      
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_build_plugin(CaloTrigger art::module
    REG_SOURCE src/CaloTrigger_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CaloCluster
      
      Offline::CalorimeterGeom
      Offline::GeometryService
      Offline::RecoDataProducts
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CaloCluster/fcl)
