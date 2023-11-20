cet_make_library(
    SOURCE
      src/BaBarMu2eField.cc
      src/DetStrawElem.cc
      src/ParticleInfo.cc
    LIBRARIES PUBLIC
      
      Offline::BFieldGeom
      Offline::ConfigTools
      Offline::DataProducts
      Offline::GeometryService
      Offline::GlobalConstantsService
      Offline::TrackerGeom
      BTrk_BaBar
      BTrk_BbrGeom
      BTrk_BField
      BTrk_DetectorModel
      BTrk_TrkBase
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
