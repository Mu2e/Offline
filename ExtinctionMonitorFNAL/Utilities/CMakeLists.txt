cet_make_library(
    SOURCE
      src/EMFRandomizationParticleDefs.cc
      src/EMFRandomizationSourceDefs.cc
      src/getCharge.cc
    LIBRARIES PUBLIC
      
      Offline::DataProducts
      Offline::GlobalConstantsService
      Offline::MCDataProducts
      Offline::ProtonBeamDumpGeom
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
