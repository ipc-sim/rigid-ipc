cmake_minimum_required(VERSION 3.1)
################################################################################


add_executable(ncp_eval ncp_eval.cpp)
target_include_directories(ncp_eval PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ncp_eval PUBLIC ${PROJECT_NAME} warnings::all)


add_executable(ncp_displ_eval ncp_displ_eval.cpp)
target_include_directories(ncp_displ_eval PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(ncp_displ_eval PUBLIC ${PROJECT_NAME} warnings::all)

if (CODE_COVERAGE)
    target_code_coverage(ncp_eval AUTO ALL)
    target_code_coverage(ncp_displ_eval AUTO ALL)
endif()
