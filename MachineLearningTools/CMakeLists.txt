cet_make_library(
    SOURCE
      src/ScoreBasedDiffusionModel.cc
    LIBRARIES PUBLIC
      CLHEP::CLHEP
      messagefacility::messagefacility
      cetlib_except::cetlib_except
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
