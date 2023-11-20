cet_make_library(
    SOURCE
      src/CosmicCORSIKA.cc
      src/ExtMonFNALMARSUtils.cc
      src/STMTestBeamFileNameDecoder.cc
      src/STMTestBeamHeaders.cc
    LIBRARIES PUBLIC
      
      Offline::CalorimeterGeom
      Offline::DataProducts
      Offline::ExtinctionMonitorFNAL_Geometry
      Offline::GeneralUtilities
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::Mu2eUtilities
      Offline::RecoDataProducts
      Offline::SeedService
)
cet_build_plugin(FromCorsikaBinary art::source
    REG_SOURCE src/FromCorsikaBinary_source.cc
    LIBRARIES REG
      Offline::Sources
      Offline::SeedService
)
cet_build_plugin(FromEMFMARSFileWeighted art::source
    REG_SOURCE src/FromEMFMARSFileWeighted_source.cc
    LIBRARIES REG
      Offline::Sources
)
cet_build_plugin(FromExtMonFNALMARSFile art::source
    REG_SOURCE src/FromExtMonFNALMARSFile_source.cc
    LIBRARIES REG
      Offline::Sources
)
cet_build_plugin(FromSTMTestBeamData art::source
    REG_SOURCE src/FromSTMTestBeamData_source.cc
    LIBRARIES REG
      Offline::Sources
)
cet_build_plugin(FromTrackerPrototypeData art::source
    REG_SOURCE src/FromTrackerPrototypeData_source.cc
    LIBRARIES REG
      Offline::Sources
      ROOT::Tree
)
cet_build_plugin(PBISequence art::source
    REG_SOURCE src/PBISequence_source.cc
    LIBRARIES REG
      Offline::Sources
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
