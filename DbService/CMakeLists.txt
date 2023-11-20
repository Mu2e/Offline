cet_make_library(
    SOURCE
      src/DbEngine.cc
      src/DbIdList.cc
      src/DbReader.cc
      src/DbService.cc
      src/DbSql.cc
      src/DbTool.cc
      src/DbValTool.cc
      src/EpicsTool.cc
    LIBRARIES PUBLIC
      PostgreSQL::PostgreSQL
      CURL::libcurl
      Offline::ConfigTools
      Offline::DbTables
      Offline::GeneralUtilities
      art::Framework_Services_Registry
)

cet_make_exec(NAME dbTool
    SOURCE src/dbTool_main.cc
    LIBRARIES
      Offline::DbService
      
)

cet_make_exec(NAME epicsTool
    SOURCE src/epicsTool_main.cc
    LIBRARIES
      Boost::program_options
      Offline::DbService
      
)

cet_build_plugin(DbService art::service
    REG_SOURCE src/DbService_service.cc
    LIBRARIES REG
      Offline::DbService
      
      Offline::ConfigTools
      Offline::DbTables
)

cet_build_plugin(DbServiceTest art::module
    REG_SOURCE src/DbServiceTest_module.cc
    LIBRARIES REG
      Offline::DbService
      
      Offline::DbTables
)

install_source(SUBDIRS src)
install_headers(SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME DbService)
