cet_make_library(
    SOURCE
      src/CosmicCORSIKA.cc
      src/ExtMonFNALMARSUtils.cc
      src/STMTestBeamFileNameDecoder.cc
      src/STMTestBeamHeaders.cc
    LIBRARIES PUBLIC
      
      offline::CalorimeterGeom
      offline::DataProducts
      offline::ExtinctionMonitorFNAL_Geometry
      offline::GeneralUtilities
      offline::GeometryService
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::Mu2eUtilities
      offline::RecoDataProducts
      offline::SeedService
)
cet_build_plugin(FromCorsikaBinary art::source
    REG_SOURCE src/FromCorsikaBinary_source.cc
    LIBRARIES REG
      offline::Sources
      offline::SeedService_SeedService_service
)
cet_build_plugin(FromEMFMARSFileWeighted art::source
    REG_SOURCE src/FromEMFMARSFileWeighted_source.cc
    LIBRARIES REG
      offline::Sources
)
cet_build_plugin(FromExtMonFNALMARSFile art::source
    REG_SOURCE src/FromExtMonFNALMARSFile_source.cc
    LIBRARIES REG
      offline::Sources
)
cet_build_plugin(FromSTMTestBeamData art::source
    REG_SOURCE src/FromSTMTestBeamData_source.cc
    LIBRARIES REG
      offline::Sources
)
cet_build_plugin(FromTrackerPrototypeData art::source
    REG_SOURCE src/FromTrackerPrototypeData_source.cc
    LIBRARIES REG
      offline::Sources
      ROOT::Tree
)
cet_build_plugin(PBISequence art::source
    REG_SOURCE src/PBISequence_source.cc
    LIBRARIES REG
      offline::Sources
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
