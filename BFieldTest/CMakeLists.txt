cet_build_plugin(BFieldSymmetry art::module
    REG_SOURCE src/BFieldSymmetry_module.cc
    LIBRARIES REG
)

cet_build_plugin(BFieldTest art::module
    REG_SOURCE src/BFieldTest_module.cc
    LIBRARIES REG
)

install_source(SUBDIRS src)
