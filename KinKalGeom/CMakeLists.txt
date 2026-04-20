cet_make_library(
    SOURCE
      src/CRV.cc
      src/DetectorSolenoid.cc
      src/StoppingTarget.cc
      src/KinKalGeom.cc
      src/TestCRV.cc
      src/Tracker.cc
    LIBRARIES PUBLIC
      KinKal::Geometry
      Offline::GeneralUtilities
      Offline::GeomPrimitives
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
