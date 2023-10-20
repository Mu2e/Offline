cet_make_library(
    SOURCE
      src/BaBarMu2eField.cc
      src/DetStrawElem.cc
      src/ParticleInfo.cc
    LIBRARIES PUBLIC
      
      offline::BFieldGeom
      offline::ConfigTools
      offline::DataProducts
      offline::GeometryService
      offline::GlobalConstantsService
      offline::TrackerGeom
      BTrk_BaBar
      BTrk_BbrGeom
      BTrk_BField
      BTrk_DetectorModel
      BTrk_TrkBase
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
