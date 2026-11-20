cet_make_library(
    SOURCE
      src/ClusterAssociator.cc
      src/ClusterFinder.cc
      src/ClusterUtils.cc
      src/GnnClusterAssembler.cc
      src/GnnGraphBuilder.cc
    LIBRARIES PUBLIC
      cetlib_except::cetlib_except
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
      Offline::ProditionsService
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

# GNN clustering: graph maker, ONNX cluster maker, parity dump, and the
# standalone assembler test.
cet_build_plugin(CaloHitGraphMaker art::module
    REG_SOURCE src/CaloHitGraphMaker_module.cc
    LIBRARIES REG
      Offline::CaloCluster
      Offline::CalorimeterGeom
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeometryService
      Offline::RecoDataProducts
)

if(onnxruntime_FOUND)
  cet_build_plugin(CaloClusterMakerGNN art::module
      REG_SOURCE src/CaloClusterMakerGNN_module.cc
      LIBRARIES REG
        onnxruntime::onnxruntime
        Offline::CaloCluster
        Offline::CalorimeterGeom
        Offline::ConfigTools
        Offline::ProditionsService
        Offline::RecoDataProducts
  )
endif()

cet_build_plugin(CaloHitGraphParityDump art::module
    REG_SOURCE src/CaloHitGraphParityDump_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      ROOT::Tree
      Offline::CalorimeterGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::RecoDataProducts
)

cet_make_exec(NAME testGnnClusterAssembler
    SOURCE src/testGnnClusterAssembler_main.cc
    LIBRARIES
      Offline::CaloCluster
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/calo_cluster_net_v2_stage1.onnx      ${CURRENT_BINARY_DIR} data/calo_cluster_net_v2_stage1.onnx      COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/calo_cluster_net_v2_stage1.norm.json ${CURRENT_BINARY_DIR} data/calo_cluster_net_v2_stage1.norm.json COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/simple_edge_net_v2.onnx              ${CURRENT_BINARY_DIR} data/simple_edge_net_v2.onnx              COPYONLY)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CaloCluster/fcl)
install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/CaloCluster)
