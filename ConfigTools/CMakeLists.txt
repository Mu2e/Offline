cet_make_library(
    SOURCE
      src/checkForStale.cc
      src/requireUniqueKey.cc
      src/SimpleConfig.cc
      src/SimpleConfigRecord.cc
    LIBRARIES PUBLIC
      messagefacility::MF_MessageLogger
      Offline::GeneralUtilities
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
