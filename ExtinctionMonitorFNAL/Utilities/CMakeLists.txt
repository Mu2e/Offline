cet_make_library(
    SOURCE
      src/EMFRandomizationParticleDefs.cc
      src/EMFRandomizationSourceDefs.cc
      src/getCharge.cc
    LIBRARIES PUBLIC
      
      offline::DataProducts
      offline::GlobalConstantsService
      offline::MCDataProducts
      offline::ProtonBeamDumpGeom
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
