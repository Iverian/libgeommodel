enable_testing()
find_package(GTest MODULE REQUIRED)

file(GLOB TEST_SOURCES "src/test_*.cpp")

# Make all tests in sigle executable
add_executable(all_tests "src/all_tests.cpp" ${TEST_SOURCES})

target_compile_definitions(all_tests PRIVATE
                           _SILENCE_ALL_CXX17_DEPRECATION_WARNINGS)

set_target_properties(all_tests
                      PROPERTIES CXX_STANDARD
                                 17
                                 CXX_STANDARD_REQUIRED
                                 YES
                                 CXX_EXTENSIONS
                                 NO
                                 RUNTIME_OUTPUT_DIRECTORY_DEBUG
                                 ${ROOT_BINARY_DIRECTORY}
                                 RUNTIME_OUTPUT_DIRECTORY_RELEASE
                                 ${ROOT_BINARY_DIRECTORY})

target_include_directories(all_tests PRIVATE "../src")

target_link_libraries(all_tests
                      PRIVATE ${PROJECT_TARGET_NAME} fmt::fmt GTest::GTest)

add_test(all_tests all_tests)
#
