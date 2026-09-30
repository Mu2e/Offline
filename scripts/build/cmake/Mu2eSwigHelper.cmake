# cmake/Mu2eSwigHelper.cmake
#
# called from a Offline subdirectory
# PYTHON_INSTALL_DIR was set at top level

function(mu2e_swig TARGET_NAME INTERFACE_FILE CPP_LIB_DEPENDENCY)
    # Define the output directory for this specific sub-module
    set(GEN_DIR "${CMAKE_CURRENT_BINARY_DIR}/swig_gen")
    set(CMAKE_SWIG_OUTDIR ${GEN_DIR})

    # Set C++ mode on the .i file
    set_property(SOURCE ${INTERFACE_FILE} PROPERTY CPLUSPLUS ON)
  
    # Create the SWIG library
    swig_add_library(${TARGET_NAME}
        LANGUAGE Python
        SOURCES ${INTERFACE_FILE}
    )

    # Inherit includes and link to the requested C++ library + Python
    set_property(TARGET ${TARGET_NAME} PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES ON)
    
    target_link_libraries(${TARGET_NAME} 
        PRIVATE 
            ${CPP_LIB_DEPENDENCY}
            Python3::Module
            Python3::Python
    )

    #target_compile_options(${TARGET_NAME} PRIVATE "-std=c++20")

    # Install both the .so and the .py file
    install(TARGETS ${TARGET_NAME}
            LIBRARY DESTINATION ${PYTHON_INSTALL_DIR})

    get_target_property(_py_dir ${TARGET_NAME} SWIG_SUPPORT_FILES_DIRECTORY)
    
    # Note: This assumes the .py filename matches the target name 
    # (or you can pass the module name as a 4th argument)
    install(DIRECTORY ${_py_dir}/
            DESTINATION ${PYTHON_INSTALL_DIR}
            FILES_MATCHING PATTERN "*.py")

endfunction()
        
