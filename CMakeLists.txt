cmake_minimum_required(VERSION 2.8)
project (graphyt)

#find_package(CUDA)
#if(CUDA_FOUND)
#
#MESSAGE(Compiling CUDA)
        # Pass options to NVCC and compile for GPU/CPU code
#        set(
#            CUDA_NVCC_FLAGS
#            ${CUDA_NVCC_FLAGS};
#            -O3 -arch=sm_35 --ptxas-options=-v
#            )#
#
#
#        # Pass options to GCC
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -O3 -std=c++11 ")
#
#
#        find_package(OpenMP REQUIRED)
#        if(OPENMP_FOUND)
#          set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#          set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#          set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#        endif()
#
#        cuda_add_library(para Parameters.cpp)
#        cuda_add_library(grid Grid.cpp)
#        cuda_add_library(extrapolation Extrapolation.cpp)
 #       cuda_add_library(TimeIntegration TimeIntegration.cu)
 #       cuda_add_library(ConstitutiveModels ConstitutiveModels.cpp)
#        cuda_add_library(DamageModels DamageModels.cpp)
#        cuda_add_library(Solver Solver.cu)
#        cuda_add_library(MaterialPoints MaterialPoints.cpp)
#        cuda_add_library(Nodes Nodes.cpp)
#        cuda_add_library(materials MaterialModels.cpp)
#        cuda_add_library(BC BoundaryConditions.cpp)
#        cuda_add_library(eig3 eig3.cpp)
#        cuda_add_library(Matrix Matrix.cpp)
#        cuda_add_executable(main main.cpp)
#
#        target_link_libraries(main MaterialPoints Nodes Solver materials para grid extrapolation TimeIntegration ConstitutiveModels DamageModels BC eig3 Matrix)
#        if(MSVC)
###          option(SWIG_DEBUG "Swig outputs to debug/ (default is release/)" OFF)
#          if(SWIG_DEBUG)
#            SET(CMAKE_SWIG_OUTDIR ${PROJECT_BINARY_DIR}/Debug)
#          else()
##            SET(CMAKE_SWIG_OUTDIR ${PROJECT_BINARY_DIR}/Release)
#          endif()
#        endif()
#        FIND_PACKAGE(SWIG REQUIRED)
#        INCLUDE(${SWIG_USE_FILE})
#        FIND_PACKAGE(PythonLibs REQUIRED)
#        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
#        INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
#        SET_SOURCE_FILES_PROPERTIES(graphyt.i PROPERTIES CPLUSPLUS ON)
#        SET_SOURCE_FILES_PROPERTIES(graphyt.i PROPERTIES SWIG_FLAGS "-includeall")
#        SWIG_ADD_MODULE(graphyt python graphyt.i)
#        SWIG_LINK_LIBRARIES(graphyt ${PYTHON_LIBRARIES})
#        SWIG_LINK_LIBRARIES(graphyt MaterialPoints Nodes Solver para grid materials extrapolation TimeIntegration ConstitutiveModels DamageModels BC eig3 Matrix)

        # option(DEBUG_MODE "DEBUG_MODE" OFF)
        # if(DEBUG_MODE)
        #    set_target_properties(_graphyt PROPERTIES DEBUG_POSTFIX  "_d")
        # endif()

#else()
MESSAGE(Compiling GCC)
        # Otherwise just compile the CPU version
#	set(CMAKE_C_COMPILER clang-omp CACHE STRING "C compiler" FORCE)
#	set(CMAKE_CXX_COMPILER clang-omp++ CACHE STRING "C++ compiler" FORCE)


        # Pass options to GCC
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpic -O3 -std=c++11 ")


        find_package(OpenMP REQUIRED)
        if(OPENMP_FOUND)
          set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
          set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
          set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
        endif()

        add_library(para Parameters.cpp)
        add_library(grid Grid.cpp)
        add_library(extrapolation Extrapolation.cpp)
        add_library(TimeIntegration TimeIntegration.cpp)
        add_library(ConstitutiveModels ConstitutiveModels.cpp)
        add_library(DamageModels DamageModels.cpp)
        add_library(Solver Solver.cpp)
        add_library(MaterialPoints MaterialPoints.cpp)
        add_library(Nodes Nodes.cpp)
        add_library(materials MaterialModels.cpp)
        add_library(BC BoundaryConditions.cpp)
        add_library(eig3 eig3.cpp)
        add_library(Matrix Matrix.cpp)
        add_executable(main main.cpp)

        target_link_libraries(main MaterialPoints Nodes Solver materials para grid extrapolation TimeIntegration ConstitutiveModels DamageModels BC eig3 Matrix)
        if(MSVC)
          option(SWIG_DEBUG "Swig outputs to debug/ (default is release/)" OFF)
          if(SWIG_DEBUG)
            SET(CMAKE_SWIG_OUTDIR ${PROJECT_BINARY_DIR}/Debug)
          else()
            SET(CMAKE_SWIG_OUTDIR ${PROJECT_BINARY_DIR}/Release)
          endif()
        endif()
        FIND_PACKAGE(SWIG REQUIRED)
        INCLUDE(${SWIG_USE_FILE})
        FIND_PACKAGE(PythonLibs REQUIRED)
        INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
        INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
        SET_SOURCE_FILES_PROPERTIES(graphyt.i PROPERTIES CPLUSPLUS ON)
        SET_SOURCE_FILES_PROPERTIES(graphyt.i PROPERTIES SWIG_FLAGS "-includeall")
        SWIG_ADD_MODULE(graphyt python graphyt.i)
        SWIG_LINK_LIBRARIES(graphyt ${PYTHON_LIBRARIES})
        SWIG_LINK_LIBRARIES(graphyt MaterialPoints Nodes Solver para grid materials extrapolation TimeIntegration ConstitutiveModels DamageModels BC eig3 Matrix)

        # option(DEBUG_MODE "DEBUG_MODE" OFF)
        # if(DEBUG_MODE)
        #    set_target_properties(_graphyt PROPERTIES DEBUG_POSTFIX  "_d")
        # endif()
#endif()

#CPACK STUFF
MESSAGE(creating installers)

set(CPACK_PACKAGE_VERSION "0.0.1")
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "1")
set(CPACK_PACKAGE_VENDOR "MIT Geonumerics")

set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_SOURCE_DIR}/LICENSE.txt")

option(USE_WIX_TOOLSET "Use WIX" true)

if (WIN32)
    if (USE_WIX_TOOLSET)
        set(CPACK_GENERATOR "WIX") # this need WiX Tooset installed and a path to candle.exe
        set(CPACK_WIX_TEMPLATE "${CMAKE_CURRENT_SOURCE_DIR}/install/windows/graphyt.template.in")
        SET(CPACK_WIX_PRODUCT_GUID F8BCDDD1-00AC-4A16-A2FE-F4E2250D2330)
        SET(CPACK_WIX_UPGRADE_GUID F8BCDDD1-00AC-4A16-A2FE-F4E2250D2330)
        # SET(CPACK_WIX_PRODUCT_ICON ../rsc/app.ico)
        SET(CPACK_WIX_UI_BANNER  "${CMAKE_CURRENT_SOURCE_DIR}/install/windows/banner_493x58.bmp")
        SET(CPACK_WIX_UI_DIALOG  "${CMAKE_CURRENT_SOURCE_DIR}/install/windows/dialog_493x312.bmp")
    else ()
        set(CPACK_GENERATOR "NSIS") # this needs NSIS installed, and available
    endif ()
    set(CPACK_SOURCE_GENERATOR "ZIP")
    set(INSTALL_DEST .)
elseif ( ${CMAKE_SYSTEM_NAME} MATCHES "Darwin" )
    set(CPACK_GENERATOR "PackageMake")
    set(INSTALL_DEST local/graphyt)

else ()
  set(CPACK_RPM_EXCLUDE_FROM_AUTO_FILELIST_ADDITION "/usr/local")
  set(INSTALL_DEST local/graphyt)
  install(TARGETS _graphyt DESTINATION ${INSTALL_DEST})
  install(FILES ${CMAKE_SOURCE_DIR}/install/unix/set_graphyt_env.sh DESTINATION ${INSTALL_DEST})
    
  option(BUILD_RPM "Builds an RPM instead of a DEB" OFF)
  if(BUILD_RPM)
    set(CPACK_RPM_PACKAGE_AUTOREQ OFF) # should probably use CPACK_RPM_PACKAGE_CONFLICTS instead
    set(CPACK_GENERATOR "RPM")
  else()
    set(CPACK_GENERATOR "DEB")
    set(CPACK_SOURCE_GENERATOR "TGZ")
    set(CPACK_DEBIAN_PACKAGE_MAINTAINER "Sam Raymond") #required
  endif()
endif ()

# INSTALL(TARGETS ${SWIG_MODULE_swig_test_REAL_NAME} DESTINATION test)
INSTALL(
  TARGETS _graphyt   
  LIBRARY DESTINATION lib  # <-- Add this line
  DESTINATION ${INSTALL_DEST})
INSTALL(FILES ${CMAKE_BINARY_DIR}/graphyt.py DESTINATION ${INSTALL_DEST})
INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/test/ DESTINATION ${INSTALL_DEST}/test FILES_MATCHING PATTERN "*.py")
INSTALL(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/Example_Scripts/ DESTINATION ${INSTALL_DEST}/Examples FILES_MATCHING PATTERN "*.py")
include(CPack)
