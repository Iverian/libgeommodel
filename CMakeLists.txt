cmake_minimum_required(VERSION 3.12)
project(libgeommodel VERSION 2.0.0 LANGUAGES CXX)

option(ENABLE_TESTS "Enable unit tests" YES)

if(EXTRA_CXX_FLAGS)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")
  message(STATUS "Added extra c++ compiler flags ${EXTRA_CXX_FLAGS}")
endif()

# Include additional CMake packages
include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
#

list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_SOURCE_DIR}/cmake")

# Configure dependencies
find_package(fmt CONFIG REQUIRED)
#

# Find source files
file(GLOB_RECURSE
  GEOMMODEL_SOURCES
  CONFIGURE_DEPENDS
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/*.hpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/include/*.hpp"
)
#

# Create target and set properties
add_library(geommodel SHARED)

target_sources(geommodel
  PRIVATE
    ${GEOMMODEL_SOURCES}
)
target_include_directories(geommodel
  PUBLIC
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_FULL_INCLUDEDIR}>"
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR};${CMAKE_CURRENT_SOURCE_DIR}/src;${CMAKE_CURRENT_SOURCE_DIR}/include>"
)
target_compile_features(geommodel
  PUBLIC
    cxx_std_17
)
target_compile_definitions(geommodel
  PUBLIC
    $<BUILD_INTERFACE:_SILENCE_ALL_CXX17_DEPRECATION_WARNINGS>
)
target_link_libraries(geommodel
  PRIVATE
    fmt::fmt
)
set_target_properties(geommodel
  PROPERTIES
    DEBUG_POSTFIX d
    EXPORT_NAME GeomModel
    CXX_EXTENSIONS NO
    WINDOWS_EXPORT_ALL_SYMBOLS ${ENABLE_TESTS}
)

include(GenerateExportHeader)
generate_export_header(geommodel
  BASE_NAME GM
  EXPORT_FILE_NAME exports.hpp
  NO_DEPRECATED_MACRO_NAME GM_OMIT_DEPRECATED
  INCLUDE_GUARD_NAME GEOMMODEL_INCLUDE_GM_EXPORTS_H_
  DEFINE_NO_DEPRECATED
)
if (ENABLE_TESTS)
  target_compile_definitions(geommodel
    PUBLIC
      $<BUILD_INTERFACE:GM_STATIC_DEFINE>
  )
endif()
#

# Install package
set(INSTALL_CONFIGDIR ${CMAKE_INSTALL_LIBDIR}/cmake/geom_model)

install(
  TARGETS geommodel
  EXPORT geommodel_targets
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}"
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  ARCHIVE DESTINATION "${CMAKE_INSTALL_LIBDIR}"
)
install(
  DIRECTORY "include/"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)
install(
  EXPORT geommodel_targets
  FILE GeomModelTargets.cmake
  NAMESPACE GeomModel::
  DESTINATION "${INSTALL_CONFIGDIR}"
)
if (WIN32)
  install(
    FILES "$<TARGET_PDB_FILE:geommodel>"
    DESTINATION "${CMAKE_INSTALL_BINDIR}"
    OPTIONAL
  )
endif()
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/exports.hpp"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/gm"
)

write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/GeomModelConfigVersion.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)
configure_package_config_file(
  "${CMAKE_CURRENT_LIST_DIR}/cmake/GeomModelConfig.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/GeomModelConfig.cmake"
  INSTALL_DESTINATION "${INSTALL_CONFIGDIR}"
)

# Install the config, configversion and custom find modules
install(
  FILES "${CMAKE_CURRENT_BINARY_DIR}/GeomModelConfig.cmake"
  "${CMAKE_CURRENT_BINARY_DIR}/GeomModelConfigVersion.cmake"
  DESTINATION "${INSTALL_CONFIGDIR}"
)
#

export(
  EXPORT geommodel_targets
  FILE "${CMAKE_CURRENT_BINARY_DIR}/GeomModelTargets.cmake"
  NAMESPACE GeomModel::
)

# Register package in user's package registry
export(PACKAGE GeomModel)
#

# Add tests
if (ENABLE_TESTS)
  add_subdirectory(tests)
endif()
#
