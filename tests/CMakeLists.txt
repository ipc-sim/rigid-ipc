# TODO: Update this code to not require installed Catch2
# https://github.com/catchorg/Catch2/blob/master/docs/cmake-integration.md
cmake_minimum_required(VERSION 3.1)

find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(Catch2 REQUIRED) # Find Catch2 for unit testing

# List all of the test source code
set(test_sources
	main.cpp
	test_ccd.cpp
	../src/collision_detection.cpp
	../src/degenerate_edge_error.cpp
	../src/not_implemented_error.cpp
	../src/impact.cpp
)

add_executable(unit_tests ${test_sources}) # Create an executable with all tests

foreach(source IN ITEMS ${test_sources})
	source_group("tests" FILES "${source}")
endforeach()

target_link_libraries(unit_tests PUBLIC Catch2::Catch2 Eigen3::Eigen)
include_directories(unit_tests ../include ../src)

# Register tests
# add_library(catch INTERFACE)
# target_include_directories(catch SYSTEM INTERFACE ${THIRD_PARTY_DIR}/catch/include)
set(PARSE_CATCH_TESTS_ADD_TO_CONFIGURE_DEPENDS ON)
include(Catch)
catch_discover_tests(unit_tests)
