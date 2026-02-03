
include(FindPackageHandleStandardArgs)
find_path(GMP_INCLUDE_DIR NAMES gmp.h HINTS ENV GMP_INC_DIR ENV GMP_DIR $ENV{GMP_DIR}/include PATH_SUFFIXES include)
find_library(GMP_LIBRARY NAMES gmp libgmp-10 gmp-10 mpir gmpxx
    HINTS ENV GMP_LIB_DIR
          ENV GMP_DIR
          $ENV{GMP_DIR}/lib
    PATH_SUFFIXES lib
    DOC "Path to the Release GMP library"
    )
find_package_handle_standard_args(GMP DEFAULT_MSG GMP_LIBRARY GMP_INCLUDE_DIR)

if(${GMP_FOUND} AND NOT TARGET GMP::GMP)
    add_library(GMP::GMP UNKNOWN IMPORTED)
    set_target_properties(GMP::GMP PROPERTIES
        IMPORTED_LINK_INTERFACE_LANGUAGES "C"
        IMPORTED_LOCATION "${GMP_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${GMP_INCLUDE_DIR}"
    )
endif()

if (NOT ${GMP_FOUND}) 
message("NOT building Blinding library because gmp-devel is not installed on this system!")
else()

include_directories(${GMP_INCLUDE_DIR})
cet_make_library(
    SOURCE
      src/BigNumber.cc
      src/GMPRoutines.cc
      src/KalSeedPrescaleTool.cc
      src/UnivariateLookupKalSeedPrescaleTool.cc
    LIBRARIES PUBLIC
      Offline::RecoDataProducts
      GMP::GMP 
)

cet_build_plugin(BigNumberChecker art::module
    REG_SOURCE src/BigNumberChecker_module.cc
    LIBRARIES REG
      Offline::Blinding
)

cet_build_plugin(BigNumberProducer art::module
    REG_SOURCE src/BigNumberProducer_module.cc
    LIBRARIES REG
      Offline::Blinding
)

cet_build_plugin(DisplaceDigiTimes art::module
    REG_SOURCE src/DisplaceDigiTimes_module.cc
    LIBRARIES REG
      Offline::Blinding
      Offline::SeedService
      Offline::Mu2eUtilities
)

cet_build_plugin(GoldwasserMicaliDecrypter art::module
    REG_SOURCE src/GoldwasserMicaliDecrypter_module.cc
    LIBRARIES REG
      Offline::Blinding
)

cet_build_plugin(GoldwasserMicaliEncrypter art::module
    REG_SOURCE src/GoldwasserMicaliEncrypter_module.cc
    LIBRARIES REG
      Offline::Blinding
      Offline::SeedService
)

cet_build_plugin(KalSeedFunctionalPrescale art::module
    REG_SOURCE src/KalSeedFunctionalPrescale_module.cc
    LIBRARIES REG
      Offline::Blinding
      Offline::SeedService
)

cet_build_plugin(MergeDigis art::module
    REG_SOURCE src/MergeDigis_module.cc
    LIBRARIES REG
      Offline::Blinding
      Offline::TrackerMC
      Offline::CaloMC
)

cet_build_plugin(TrackDigiExtractor art::module
    REG_SOURCE src/TrackDigiExtractor_module.cc
    LIBRARIES REG
      Offline::Blinding
)

cet_build_plugin(MomentumLookupKalSeedPrescaleTool art::tool
    REG_SOURCE src/MomentumLookupKalSeedPrescaleTool_tool.cc
    LIBRARIES REG
        Offline::Blinding
)

cet_build_plugin(QuasiImpactParameterLookupKalSeedPrescaleTool art::tool
    REG_SOURCE src/QuasiImpactParameterLookupKalSeedPrescaleTool_tool.cc
    LIBRARIES REG
        Offline::Blinding
)

cet_build_plugin(TrackerExitAzimuthLookupKalSeedPrescaleTool art::tool
    REG_SOURCE src/TrackerExitAzimuthLookupKalSeedPrescaleTool_tool.cc
    LIBRARIES REG
        Offline::Blinding
)

art_dictionary( NO_CHECK_CLASS_VERSION # For some reason this segfaults
    CLASSES_DEF_XML ${CMAKE_CURRENT_SOURCE_DIR}/src/classes_def.xml
    CLASSES_H ${CMAKE_CURRENT_SOURCE_DIR}/src/classes.h
     DICTIONARY_LIBRARIES
      Offline::Blinding
)

endif()

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
