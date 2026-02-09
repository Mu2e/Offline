cet_build_plugin(CaloClusterCheck art::module
    REG_SOURCE src/CaloClusterCheck_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CaloExample art::module
    REG_SOURCE src/CaloExample_module.cc
    LIBRARIES REG      
      art_root_io::TFileService_service
      Offline::CaloCluster
      Offline::CalorimeterGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CaloMCInspector art::module
    REG_SOURCE src/CaloMCInspector_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CaloNeutron art::module
    REG_SOURCE src/CaloNeutron_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
)

cet_build_plugin(CaloNNDiag art::module
    REG_SOURCE src/CaloNNDiag_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CaloCluster
      Offline::CalorimeterGeom
      Offline::ConditionsService
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

install_source(SUBDIRS src)
