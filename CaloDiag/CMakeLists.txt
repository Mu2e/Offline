cet_build_plugin(CaloClusterCheck art::module
    REG_SOURCE src/CaloClusterCheck_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::CalorimeterGeom
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(CaloExample art::module
    REG_SOURCE src/CaloExample_module.cc
    LIBRARIES REG      
      art_root_io::TFileService_service
      offline::CaloCluster
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(CaloMCInspector art::module
    REG_SOURCE src/CaloMCInspector_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::MCDataProducts
      offline::RecoDataProducts
)

cet_build_plugin(CaloNeutron art::module
    REG_SOURCE src/CaloNeutron_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
)

cet_build_plugin(CaloNNDiag art::module
    REG_SOURCE src/CaloNNDiag_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      offline::CaloCluster
      offline::CalorimeterGeom
      offline::ConditionsService
      offline::DataProducts
      offline::GeometryService
      offline::MCDataProducts
      offline::RecoDataProducts
)

install_source(SUBDIRS src)
