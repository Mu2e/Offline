cet_make_library(
    SOURCE
      src/DbEngine.cc
      src/DbIdList.cc
      src/DbReader.cc
      src/DbSql.cc
      src/DbTool.cc
      src/DbValTool.cc
      src/EpicsTool.cc
    LIBRARIES PUBLIC
      PostgreSQL::PostgreSQL
      CURL::libcurl
      offline::ConfigTools
      offline::DbTables
      offline::GeneralUtilities
)

cet_make_exec(NAME dbTool
    SOURCE src/dbTool_main.cc 
    LIBRARIES offline::DbService
)

cet_make_exec(NAME epicsTool 
    SOURCE src/epicsTool_main.cc 
    LIBRARIES 
      offline::DbService
      Boost::program_options
)

cet_build_plugin(DbService art::service
    REG_SOURCE src/DbService_service.cc
    LIBRARIES REG
      offline::DbService
)

cet_build_plugin(DbServiceTest art::module
    REG_SOURCE src/DbServiceTest_module.cc
    LIBRARIES REG
      offline::DbService
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME DbService)
