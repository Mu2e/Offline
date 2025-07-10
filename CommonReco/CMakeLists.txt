cet_build_plugin(CosmicTimeOffset art::module
  REG_SOURCE src/SelectReco_module.cc
  LIBRARIES REG
  Offline::ConditionsService
  Offline::DataProducts
  Offline::GeometryService
  Offline::ProditionsService
  Offline::RecoDataProducts
  Offline::TrackerConditions
  Offline::TrackerGeom
  )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)

install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CommonReco/fcl)
