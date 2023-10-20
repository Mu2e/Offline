cet_make_library(
    SOURCE
      src/DetectorSolenoid.cc
      src/StoppingTarget.cc
      src/SurfaceId.cc
      src/SurfaceMap.cc
      src/TestCRV.cc
      src/Tracker.cc
    LIBRARIES PUBLIC
      KinKal_Geometry
      offline::GeneralUtilities
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
