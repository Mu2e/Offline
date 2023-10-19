cet_make_library(
    SOURCE
src/DbEngine.cc
src/DbIdList.cc
src/DbReader.cc
src/DbService_service.cc
src/DbSql.cc
src/DbTool.cc
src/dbTool_main.cc
src/DbValTool.cc
src/EpicsTool.cc
src/epicsTool_main.cc
LIBRARIES PUBLIC
cetlib::cetlib
)

cet_build_plugin(DbServiceTest art::module REG_SOURCE src/DbServiceTest_module.cc)


install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME DbService)
