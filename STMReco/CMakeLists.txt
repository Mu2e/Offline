cet_build_plugin(MakeSTMHits art::module
    REG_SOURCE src/MakeSTMHits_module.cc
    LIBRARIES REG
      Offline::GlobalConstantsService
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::STMConditions
)

cet_build_plugin(PlotSTMEnergySpectrum art::module
    REG_SOURCE src/PlotSTMEnergySpectrum_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(PlotSTMMWDSpectrum art::module
    REG_SOURCE src/PlotSTMMWDSpectrum_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(PlotSTMWaveformDigis art::module
    REG_SOURCE src/PlotSTMWaveformDigis_module.cc
    LIBRARIES REG
      Offline::DataProducts
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::STMConditions
)

cet_build_plugin(STMMovingWindowDeconvolution art::module
    REG_SOURCE src/STMMovingWindowDeconvolution_module.cc
    LIBRARIES REG
      Offline::GlobalConstantsService
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::STMConditions
)

cet_build_plugin(STMZeroSuppression art::module
    REG_SOURCE src/STMZeroSuppression_module.cc
    LIBRARIES REG
      Offline::Mu2eUtilities
      Offline::ProditionsService
      Offline::RecoDataProducts
      Offline::STMConditions
)

cet_build_plugin(STMAnalyzeDigis art::module
    REG_SOURCE src/STMAnalyzeDigis_module.cc
    LIBRARIES REG
    art_root_io::TFileService_service
    artdaq_core_mu2e::artdaq-core-mu2e_Data
    otsdaq_mu2e::otsdaq-mu2e_ArtModules
    Offline::Mu2eUtilities
    Offline::ProditionsService
    Offline::DataProducts
    Offline::RecoDataProducts
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/makeSTMHits.fcl   ${CURRENT_BINARY_DIR} fcl/makeSTMHits.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/makeSTMHits_testbeam.fcl   ${CURRENT_BINARY_DIR} fcl/makeSTMHits_testbeam.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/mwd.fcl   ${CURRENT_BINARY_DIR} fcl/mwd.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/mwd_testbeam.fcl   ${CURRENT_BINARY_DIR} fcl/mwd_testbeam.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/plotSTMEnergySpectrum.fcl   ${CURRENT_BINARY_DIR} fcl/plotSTMEnergySpectrum.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/plotSTMMWDSpectrum.fcl   ${CURRENT_BINARY_DIR} fcl/plotSTMMWDSpectrum.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/plotSTMWaveformDigis.fcl   ${CURRENT_BINARY_DIR} fcl/plotSTMWaveformDigis.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog_testbeam.fcl   ${CURRENT_BINARY_DIR} fcl/prolog_testbeam.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/zeroSuppression.fcl   ${CURRENT_BINARY_DIR} fcl/zeroSuppression.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/zeroSuppression_testbeam.fcl   ${CURRENT_BINARY_DIR} fcl/zeroSuppression_testbeam.fcl   COPYONLY)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/analyzeDigis.fcl   ${CURRENT_BINARY_DIR} fcl/analyzeDigis.fcl   COPYONLY)


install_source(SUBDIRS src)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/STMReco/fcl)
