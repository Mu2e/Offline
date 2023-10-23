cet_make_library(
    SOURCE
      src/ClusterAssociator.cc
      src/ClusterFinder.cc
      src/ClusterUtils.cc
    LIBRARIES PUBLIC
      
      offline::CalorimeterGeom
      offline::GeometryService
      offline::RecoDataProducts
)

cet_build_plugin(CaloClusterFast art::module
    REG_SOURCE src/CaloClusterFast_module.cc
    LIBRARIES REG
      offline::CaloCluster
      
      offline::CalorimeterGeom
      offline::GeometryService
      offline::RecoDataProducts
)

cet_build_plugin(CaloClusterMaker art::module
    REG_SOURCE src/CaloClusterMaker_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::CaloCluster
      
      offline::CalorimeterGeom
      offline::GeometryService
      offline::RecoDataProducts
)

cet_build_plugin(CaloProtoClusterMaker art::module
    REG_SOURCE src/CaloProtoClusterMaker_module.cc
    LIBRARIES REG
      offline::CaloCluster
      
      offline::CalorimeterGeom
      offline::GeometryService
      offline::RecoDataProducts
)

cet_build_plugin(CaloTrigger art::module
    REG_SOURCE src/CaloTrigger_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::CaloCluster
      
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::GeometryService
      offline::RecoDataProducts
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME CaloCluster)
