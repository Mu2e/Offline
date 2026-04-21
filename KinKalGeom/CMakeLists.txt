cet_make_library(
    SOURCE
      src/KinKalGeom.cc
    LIBRARIES PUBLIC
      KinKal::Geometry
      KinKal::Trajectory
      KinKal::General
      Offline::GeneralUtilities
      Offline::GeomPrimitives
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
