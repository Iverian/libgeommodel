cmake_minimum_required(VERSION 3.0)
project(libgeommodel VERSION 1.0.0 LANGUAGES CXX)

if(${CMAKE_BUILD_TYPE} STREQUAL "Debug" OR ${CMAKE_BUILD_TYPE} STREQUAL "RelWithDebInfo")
    set(PROJECT_TARGET_NAME "geommodeld")
else()
    set(PROJECT_TARGET_NAME "geommodel")
endif()


####################################################
# Include additional cmake packages
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
####################################################

list(INSERT CMAKE_MODULE_PATH 0 ${CMAKE_SOURCE_DIR}/cmake)

####################################################
# Configure conan
if(NOT EXISTS ${CMAKE_BINARY_DIR}/conan.cmake)
    message(STATUS "Downloading conan.cmake from https://github.com/conan-io/cmake-conan")
    file(DOWNLOAD "https://raw.githubusercontent.com/conan-io/cmake-conan/master/conan.cmake"
                  "${CMAKE_BINARY_DIR}/conan.cmake")
endif()

include(${CMAKE_BINARY_DIR}/conan.cmake)

conan_cmake_run(
    REQUIRES gtest/1.8.1@bincrafters/stable
             fmt/5.2.1@bincrafters/stable
    BASIC_SETUP CMAKE_TARGETS
    BUILD missing
)

set(CMAKE_CXX_STANDARD 17)
####################################################

####################################################
# Find source files
file(GLOB MAIN_HEADERS "include/gm/*.h")

set(MAIN_SOURCES "")
file(GLOB _SRC_CHILDREN "src/*")

foreach(_CHILD ${_SRC_CHILDREN})
    if(IS_DIRECTORY ${_CHILD})
        file(GLOB _CHILD_HEADERS "${_CHILD}/*.h")
	    file(GLOB _CHILD_SOURCES "${_CHILD}/*.cpp")
        list(APPEND MAIN_SOURCES ${_CHILD_HEADERS})
        list(APPEND MAIN_SOURCES ${_CHILD_SOURCES})
	endif()
endforeach()
####################################################

####################################################
# Create target and set properties
if(WIN32)
    set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

add_library(
    ${PROJECT_TARGET_NAME}
    SHARED
        ${MAIN_HEADERS}
        ${MAIN_SOURCES}
)

target_include_directories(
    ${PROJECT_TARGET_NAME}
    PUBLIC
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_FULL_INCLUDEDIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        "src"
)

target_link_libraries(
    ${PROJECT_TARGET_NAME}
    PRIVATE
        CONAN_PKG::fmt
)
####################################################

####################################################
# Install package
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/geom_model)

install(
    TARGETS ${PROJECT_TARGET_NAME}
    EXPORT ${PROJECT_TARGET_NAME}-targets
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
)

set_target_properties(${PROJECT_TARGET_NAME} PROPERTIES EXPORT_NAME GeomModel)

install(
    DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(
    EXPORT ${PROJECT_TARGET_NAME}-targets
    FILE GeomModelTargets.cmake
    NAMESPACE GeomModel::
    DESTINATION ${INSTALL_CONFIGDIR}
)

install(
    FILES $<TARGET_PDB_FILE:${PROJECT_TARGET_NAME}>
    DESTINATION ${CMAKE_INSTALL_BINDIR}
    OPTIONAL
)

write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/GeomModelConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

configure_package_config_file(
    ${CMAKE_CURRENT_LIST_DIR}/cmake/GeomModelConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/GeomModelConfig.cmake
    INSTALL_DESTINATION ${INSTALL_CONFIGDIR}
)

#Install the config, configversion and custom find modules
install(
    FILES
        ${CMAKE_CURRENT_BINARY_DIR}/GeomModelConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/GeomModelConfigVersion.cmake
    DESTINATION ${INSTALL_CONFIGDIR}
)
####################################################

export(
    EXPORT ${PROJECT_TARGET_NAME}-targets
    FILE ${CMAKE_CURRENT_BINARY_DIR}/GeomModelTargets.cmake
    NAMESPACE GeomModel::
)

#Register package in user's package registry
export(PACKAGE GeomModel)

####################################################
# Add tests
add_subdirectory(test)
####################################################
