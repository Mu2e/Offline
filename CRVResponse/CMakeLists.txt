cet_make_library(
    SOURCE
      src/CrvMCHelper.cc
      src/DrawLookupTableHistograms.cc
      src/MakeCrvDigis.cc
      src/MakeCrvPhotons.cc
      src/MakeCrvSiPMCharges.cc
      src/MakeCrvWaveforms.cc
    LIBRARIES PUBLIC

      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CrvCoincidenceClusterMatchMC art::module
    REG_SOURCE src/CrvCoincidenceClusterMatchMC_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CRVResponse

      Offline::ConditionsService
      Offline::CosmicRayShieldGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CrvDigitizer art::module
    REG_SOURCE src/CrvDigitizer_module.cc
    LIBRARIES REG
      Offline::CRVResponse

      Offline::CosmicRayShieldGeom
      Offline::CRVConditions
      Offline::DAQConditions
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::RecoDataProducts
)

cet_build_plugin(CrvPhotonGenerator art::module
    REG_SOURCE src/CrvPhotonGenerator_module.cc
    LIBRARIES REG
      Offline::CRVResponse

      Offline::ConditionsService
      Offline::ConfigTools
      Offline::CosmicRayShieldGeom
      Offline::DAQConditions
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::SeedService
)

cet_build_plugin(CrvPlot art::module
    REG_SOURCE src/CrvPlot_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CRVResponse

      Offline::ConditionsService
      Offline::CosmicRayShieldGeom
      Offline::CRVConditions
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::RecoDataProducts
)

cet_build_plugin(CrvSiPMChargeGenerator art::module
    REG_SOURCE src/CrvSiPMChargeGenerator_module.cc
    LIBRARIES REG
      Offline::CRVResponse

      Offline::ConditionsService
      Offline::ConfigTools
      Offline::CosmicRayShieldGeom
      Offline::DAQConditions
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::SeedService
)

cet_build_plugin(CrvStepsFromStepPointMCs art::module
    REG_SOURCE src/CrvStepsFromStepPointMCs_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CRVResponse

      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::ProditionsService
)

cet_build_plugin(CRVTest art::module
    REG_SOURCE src/CRVTest_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::CRVResponse

      Offline::CosmicRayShieldGeom
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::RecoDataProducts
)

cet_build_plugin(CrvWaveformsGenerator art::module
    REG_SOURCE src/CrvWaveformsGenerator_module.cc
    LIBRARIES REG
      Offline::CRVResponse

      Offline::ConfigTools
      Offline::CosmicRayShieldGeom
      Offline::CRVConditions
      Offline::DAQConditions
      Offline::DataProducts
      Offline::GeometryService
      Offline::MCDataProducts
      Offline::ProditionsService
      Offline::SeedService
)

cet_build_plugin(CrvWidebandTest art::module
    REG_SOURCE src/CrvWidebandTest_module.cc
    LIBRARIES REG
      Offline::CRVResponse
      art_root_io::TFileService_service
)

cet_build_plugin(MakeCrvCoincidenceClusterMCAssns art::module
    REG_SOURCE src/MakeCrvCoincidenceClusterMCAssns_module.cc
    LIBRARIES REG
      Offline::CRVResponse

      Offline::MCDataProducts
)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/epilog_extracted.fcl	  ${CURRENT_BINARY_DIR} fcl/epilog_extracted.fcl	 )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/epilog_extracted_v02.fcl   ${CURRENT_BINARY_DIR} fcl/epilog_extracted_v02.fcl )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/epilog_extracted_v03.fcl   ${CURRENT_BINARY_DIR} fcl/epilog_extracted_v03.fcl )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog_v03.fcl		  ${CURRENT_BINARY_DIR} fcl/prolog_v03.fcl		 )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog_v05.fcl		  ${CURRENT_BINARY_DIR} fcl/prolog_v05.fcl		 )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog_v06.fcl		  ${CURRENT_BINARY_DIR} fcl/prolog_v06.fcl		 )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog_v07.fcl		  ${CURRENT_BINARY_DIR} fcl/prolog_v07.fcl		 )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog_v08.fcl		  ${CURRENT_BINARY_DIR} fcl/prolog_v08.fcl		 )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog_v09.fcl		  ${CURRENT_BINARY_DIR} fcl/prolog_v09.fcl		 )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog_v10.fcl		  ${CURRENT_BINARY_DIR} fcl/prolog_v10.fcl		 )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog_v11.fcl		  ${CURRENT_BINARY_DIR} fcl/prolog_v11.fcl		 )

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/singlePEWaveform.txt	  ${CURRENT_BINARY_DIR} data/singlePEWaveform.txt	 )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/singlePEWaveform_v2.txt	  ${CURRENT_BINARY_DIR} data/singlePEWaveform_v2.txt	 )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/singlePEWaveform_v3.txt   ${CURRENT_BINARY_DIR} data/singlePEWaveform_v3.txt      )

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/CRVResponse)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/CRVResponse/fcl)
