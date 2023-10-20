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
      offline::DataProducts
      offline::GeneralUtilities
      offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
