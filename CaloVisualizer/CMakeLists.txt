ROOT_GENERATE_DICTIONARY(THMu2eCaloDisk_dict
 THMu2eCaloDisk.h
 LINKDEF LinkDef.h
 MODULE THMu2eCaloDisk
)

set_property(SOURCE ${CMAKE_CURRENT_BINARY_DIR}/THMu2eCaloDisk_dict.cxx PROPERTY GENERATED TRUE)
add_library(THMu2eCaloDisk SHARED THMu2eCaloDisk.cxx ${CMAKE_CURRENT_BINARY_DIR}/THMu2eCaloDisk_dict.cxx)
target_include_directories(THMu2eCaloDisk PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(THMu2eCaloDisk ROOT::Core ROOT::Hist ROOT::RIO Offline::DataProducts)

install(TARGETS THMu2eCaloDisk)
install_headers()
install_source()
