cet_make_library(
    SOURCE
      src/TestClass.cc
    LIBRARIES PUBLIC
      
)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
