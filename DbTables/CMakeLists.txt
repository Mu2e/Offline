cet_make_library(
    SOURCE
      src/DbCache.cc
      src/DbIoV.cc
      src/DbSet.cc
      src/DbTable.cc
      src/DbTableFactory.cc
      src/DbUtil.cc
      src/DbValCache.cc
    LIBRARIES PUBLIC
      
      Offline::DataProducts
      Offline::GeneralUtilities
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
