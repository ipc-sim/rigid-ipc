################################################################################
# MinDistance Executable
################################################################################
add_executable(cli_mindistance cli_min_distance.cpp)

target_link_libraries(cli_mindistance PUBLIC ipc::rigid)

include(cli11)
target_link_libraries(cli_mindistance PUBLIC CLI11::CLI11)

set_target_properties(cli_mindistance PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

################################################################################
# CCD Check
################################################################################
add_executable(cli_ccd cli_ccd.cpp)

target_link_libraries(cli_ccd PUBLIC ipc::rigid)

include(cli11)
target_link_libraries(cli_ccd PUBLIC CLI11::CLI11)

set_target_properties(cli_ccd PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

################################################################################
# Headless Renderer
################################################################################
add_executable(render_simulation
  renderer/main.cpp
)

target_link_libraries(render_simulation PUBLIC ipc::rigid)

include(software_renderer)
target_link_libraries(render_simulation PUBLIC software_renderer::software_renderer)

include(cli11)
target_link_libraries(render_simulation PUBLIC CLI11::CLI11)

set_target_properties(render_simulation PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

################################################################################
# OBJ Sequence Generator
################################################################################
add_executable(obj_sequence obj_sequence.cpp)

target_link_libraries(obj_sequence PUBLIC ipc::rigid)

include(cli11)
target_link_libraries(obj_sequence PUBLIC CLI11::CLI11)

set_target_properties(obj_sequence PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

################################################################################
# GLTF Generator
################################################################################
add_executable(sim_to_gltf sim_to_gltf.cpp)

target_link_libraries(sim_to_gltf PUBLIC ipc::rigid)

include(cli11)
target_link_libraries(sim_to_gltf PUBLIC CLI11::CLI11)

set_target_properties(sim_to_gltf PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

################################################################################
# CCD Comparisons
################################################################################
add_executable(ccd_comparison time_ccd.cpp)

target_link_libraries(ccd_comparison PUBLIC ipc::rigid)

include(cli11)
target_link_libraries(ccd_comparison PUBLIC CLI11::CLI11)

set_target_properties(ccd_comparison PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

################################################################################
# JSON to MJCF
################################################################################
add_executable(json_to_mjcf json_to_mjcf.cpp)

target_link_libraries(json_to_mjcf PUBLIC ipc::rigid)

include(cli11)
target_link_libraries(json_to_mjcf PUBLIC CLI11::CLI11)

include(tinyxml2)
target_link_libraries(json_to_mjcf PUBLIC tinyxml2::tinyxml2)

set_target_properties(json_to_mjcf PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

################################################################################
# Generate Bullet results
################################################################################
add_executable(generate_bullet_results generate_bullet_results.cpp)

target_link_libraries(generate_bullet_results PUBLIC ipc::rigid)

include(cli11)
target_link_libraries(generate_bullet_results PUBLIC CLI11::CLI11)

set_target_properties(generate_bullet_results PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
