cet_make_library(
    SOURCE
src/DbCache.cc
src/DbIoV.cc
src/DbSet.cc
src/DbTable.cc
src/DbTableFactory.cc
src/DbUtil.cc
src/DbValCache.cc
LIBRARIES PUBLIC
cetlib::cetlib
)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
