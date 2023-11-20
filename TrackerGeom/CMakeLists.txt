cet_make_library(
    SOURCE
      src/Manifold.cc
      src/Panel.cc
      src/Plane.cc
      src/Straw.cc
      src/SupportModel.cc
      src/SupportStructure.cc
      src/Tracker.cc
    LIBRARIES PUBLIC
      
      Offline::DataProducts
      Offline::GeneralUtilities
      Offline::GeomPrimitives
      Offline::Mu2eInterfaces
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
