cet_build_plugin(BFieldSymmetry art::module
    REG_SOURCE src/BFieldSymmetry_module.cc
    LIBRARIES REG
      offline::BFieldGeom
      offline::GeometryService
      offline::SeedService
)

cet_build_plugin(BFieldTest art::module
    REG_SOURCE src/BFieldTest_module.cc
    LIBRARIES REG
      offline::BFieldGeom
      offline::GeneralUtilities
      offline::GeometryService
)

install_source(SUBDIRS src)
