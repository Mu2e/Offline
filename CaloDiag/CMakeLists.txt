cet_build_plugin(CaloClusterCheck art::module REG_SOURCE src/CaloClusterCheck_module.cc)

cet_build_plugin(CaloExample art::module REG_SOURCE src/CaloExample_module.cc)

cet_build_plugin(CaloMCInspector art::module REG_SOURCE src/CaloMCInspector_module.cc)

cet_build_plugin(CaloNeutron art::module REG_SOURCE src/CaloNeutron_module.cc)

cet_build_plugin(CaloNNDiag art::module REG_SOURCE src/CaloNNDiag_module.cc)


install_source(SUBDIRS src)
