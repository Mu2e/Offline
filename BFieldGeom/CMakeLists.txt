cet_make_library(
    SOURCE
src/BFCacheManager.cc
src/BFGridMap.cc
src/BFieldManager.cc
src/BFInterpolationStyle.cc
src/BFMapType.cc
src/BFParamMap.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
