cet_build_plugin(HelloProducer art::module REG_SOURCE src/HelloProducer_module.cc)

cet_build_plugin(HelloWorld2 art::module REG_SOURCE src/HelloWorld2_module.cc)

cet_build_plugin(HelloWorld art::module REG_SOURCE src/HelloWorld_module.cc)


install_source(SUBDIRS src)
