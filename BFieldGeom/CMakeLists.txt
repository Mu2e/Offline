cet_make_library(
    SOURCE
      src/BFCacheManager.cc
      src/BFGridMap.cc
      src/BFieldManager.cc
      src/BFInterpolationStyle.cc
      src/BFMapType.cc
      src/BFParamMap.cc
    LIBRARIES PUBLIC
      messagefacility::MF_MessageLogger
      GSL::gsl
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
