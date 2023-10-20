cet_build_plugin(ModifyTrackSPM art::module
    REG_SOURCE src/ModifyTrackSPM_module.cc
    LIBRARIES REG
)

cet_build_plugin(NullProducer art::module
    REG_SOURCE src/NullProducer_module.cc
    LIBRARIES REG
)

cet_build_plugin(ReadTrackSPM art::module
    REG_SOURCE src/ReadTrackSPM_module.cc
    LIBRARIES REG
)

install_source(SUBDIRS src)
