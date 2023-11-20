cet_build_plugin(BFieldSymmetry art::module
    REG_SOURCE src/BFieldSymmetry_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::BFieldGeom
      Offline::GeometryService
      Offline::SeedService
      ROOT::Hist
      ROOT::Tree
)

cet_build_plugin(BFieldTest art::module
    REG_SOURCE src/BFieldTest_module.cc
    LIBRARIES REG
      Offline::BFieldGeom
      Offline::GeneralUtilities
      Offline::GeometryService
)

install_source(SUBDIRS src)
