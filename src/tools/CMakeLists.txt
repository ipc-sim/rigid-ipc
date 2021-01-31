cmake_minimum_required(VERSION 3.8)


################################################################################
# MinDistance Executable
################################################################################
add_executable(cli_mindistance
  cli_min_distance.cpp
)

target_link_libraries(cli_mindistance PUBLIC ${PROJECT_NAME})
target_link_libraries(cli_mindistance PUBLIC CLI11::CLI11)

set_target_properties(cli_mindistance PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

################################################################################
# CCD Check
################################################################################
add_executable(cli_ccd
  cli_ccd.cpp
)

target_link_libraries(cli_ccd PUBLIC ${PROJECT_NAME})
target_link_libraries(cli_ccd PUBLIC CLI11::CLI11)

set_target_properties(cli_ccd PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

################################################################################
# Headless Renderer
################################################################################
add_library(software_renderer
  renderer/camera.cpp
  renderer/raster.cpp
  renderer/render_mesh.cpp
  renderer/write_png.cpp
)

target_link_libraries(software_renderer PUBLIC igl::core igl_stb_image)
target_link_libraries(software_renderer PUBLIC nlohmann::json)
target_link_libraries(software_renderer PUBLIC TBB::tbb)

add_executable(render_simulation
  renderer/main.cpp
)

target_link_libraries(render_simulation PUBLIC software_renderer)
target_link_libraries(render_simulation PUBLIC ${PROJECT_NAME})
target_link_libraries(render_simulation PUBLIC CLI11::CLI11)

set_target_properties(render_simulation PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

################################################################################
# OBJ Sequence Generator
################################################################################
add_executable(obj_sequence
  obj_sequence.cpp
)

target_link_libraries(obj_sequence PUBLIC ${PROJECT_NAME})
target_link_libraries(obj_sequence PUBLIC CLI11::CLI11)

set_target_properties(obj_sequence PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

################################################################################
# GLTF Generator
################################################################################
add_executable(sim_to_gltf
  sim_to_gltf.cpp
)

target_link_libraries(sim_to_gltf PUBLIC ${PROJECT_NAME})
target_link_libraries(sim_to_gltf PUBLIC CLI11::CLI11)

set_target_properties(sim_to_gltf PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")

################################################################################
# CCD Comparisons
################################################################################
add_executable(ccd_comparison
  time_ccd.cpp
)

target_link_libraries(ccd_comparison PUBLIC ${PROJECT_NAME})
target_link_libraries(ccd_comparison PUBLIC CLI11::CLI11)

set_target_properties(ccd_comparison PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tools")
