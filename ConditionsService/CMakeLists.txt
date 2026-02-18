cet_make_library(
    SOURCE
      src/ConditionsService.cc
    LIBRARIES PUBLIC
      
      Offline::ConfigTools
      Offline::Mu2eInterfaces
      Offline::GeometryService
      art::Framework_Services_Registry
)

cet_build_plugin(ConditionsService art::service
    REG_SOURCE src/ConditionsService_service.cc
    LIBRARIES REG
      Offline::ConditionsService
      
)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/conditions_01.txt		      ${CURRENT_BINARY_DIR} data/conditions_01.txt		     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/pid_ele_dt_v4_2_4.tab		      ${CURRENT_BINARY_DIR} data/pid_ele_dt_v4_2_4.tab		     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/pid_ele_ep_vs_path_v4_2_4.tab	      ${CURRENT_BINARY_DIR} data/pid_ele_ep_vs_path_v4_2_4.tab	     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/pid_ele_xs_v4_2_4.tab		      ${CURRENT_BINARY_DIR} data/pid_ele_xs_v4_2_4.tab		     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/pid_muo_dt_v4_2_4.tab		      ${CURRENT_BINARY_DIR} data/pid_muo_dt_v4_2_4.tab		     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/pid_muo_ep_vs_path_v4_2_4.tab	      ${CURRENT_BINARY_DIR} data/pid_muo_ep_vs_path_v4_2_4.tab	     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/pid_muo_xs_v4_2_4.tab		      ${CURRENT_BINARY_DIR} data/pid_muo_xs_v4_2_4.tab		     )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_9/pid_ele_dedx.rtbl       ${CURRENT_BINARY_DIR} data/v5_7_9/pid_ele_dedx.rtbl      )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_9/pid_ele_dt.tbl	       ${CURRENT_BINARY_DIR} data/v5_7_9/pid_ele_dt.tbl	        )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_9/pid_ele_ep_vs_path.tbl  ${CURRENT_BINARY_DIR} data/v5_7_9/pid_ele_ep_vs_path.tbl )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_9/pid_ele_xdrds.tbl       ${CURRENT_BINARY_DIR} data/v5_7_9/pid_ele_xdrds.tbl      )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_9/pid_muo_dedx.rtbl       ${CURRENT_BINARY_DIR} data/v5_7_9/pid_muo_dedx.rtbl      )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_9/pid_muo_dt.tbl	       ${CURRENT_BINARY_DIR} data/v5_7_9/pid_muo_dt.tbl	        )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_9/pid_muo_ep_vs_path.tbl  ${CURRENT_BINARY_DIR} data/v5_7_9/pid_muo_ep_vs_path.tbl )
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/v5_7_9/pid_muo_xdrds.tbl       ${CURRENT_BINARY_DIR} data/v5_7_9/pid_muo_xdrds.tbl      )

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/ConditionsService)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
