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
      src/GrlHeader.cc
      src/GrlList.cc
      src/GrlTool.cc
      src/OpenTool.cc
      src/RunInfo.cc
      src/RunTool.cc
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

cet_make_exec(NAME runTool
    SOURCE src/runTool_main.cc
    LIBRARIES
      Boost::program_options
      Offline::DbService      
)

cet_make_exec(NAME grlTool
    SOURCE src/grlTool_main.cc
    LIBRARIES
      Boost::program_options
      Offline::DbService      
)

cet_make_exec(NAME openTool
    SOURCE src/openTool_main.cc
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

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl ${CURRENT_BINARY_DIR} fcl/prolog.fcl)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/connections.txt ${CURRENT_BINARY_DIR} data/connections.txt)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/create.sql ${CURRENT_BINARY_DIR} data/create.sql)

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/DbService)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/DbService/fcl)
