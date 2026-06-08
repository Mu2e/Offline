cet_make_library(
    SOURCE
      src/KinKalGeom.cc
      src/KKMaterial.cc
      src/KKStrawMaterial.cc
    LIBRARIES PUBLIC
      Offline::ConfigTools
      KinKal::MatEnv
      KinKal::Geometry
      KinKal::Trajectory
      KinKal::General
      Offline::GeneralUtilities
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/KinKalGeom/fcl)
