enable_testing()


file(GLOB TEST_SOURCES "src/test_*.cpp")

####################################################
# Make all tests in sigle executable
add_executable(all_tests
    "src/all_tests.cpp"
    ${TEST_SOURCES}
)

target_include_directories(all_tests
    PRIVATE
        "../src"
        ${CONAN_INCLUDE_DIRS_GTEST}
)

target_link_libraries(all_tests
    PRIVATE
        ${PROJECT_TARGET_NAME}
        CONAN_PKG::gtest
        CONAN_PKG::fmt
)

add_test(all_tests all_tests)
####################################################
