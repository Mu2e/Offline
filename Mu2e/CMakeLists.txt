cet_make_exec(NAME mu2e
    SOURCE src/mu2e_main.cc
    LIBRARIES
      art::Framework_Core
      art::Framework_Art
)

install_source(SUBDIRS src)
