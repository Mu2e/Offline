cet_make_library(
    SOURCE
      src/PIDLogL1D.cc
      src/PIDLogLEp.cc
      src/PIDLogLRatio.cc
      src/PIDUtilities.cc
    LIBRARIES PUBLIC
      ROOT::Hist
      Offline::ConfigTools
      Offline::GeneralUtilities
)

cet_build_plugin(ParticleIDRead art::module
    REG_SOURCE src/ParticleIDRead_module.cc
    LIBRARIES REG
      art_root_io::TFileService_service
      Offline::ParticleID
      
      Offline::RecoDataProducts
)


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/fcl/prolog.fcl   ${CURRENT_BINARY_DIR} fcl/prolog.fcl   COPYONLY)

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

install(DIRECTORY data DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/Offline/ParticleID)

install_source(SUBDIRS src)
install_headers(USE_PROJECT_NAME SUBDIRS inc)
install_fhicl(SUBDIRS fcl SUBDIRNAME Offline/ParticleID/fcl)
