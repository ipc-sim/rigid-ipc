cmake_minimum_required(VERSION 3.1)
################################################################################

add_executable(Box2D-comparison
    src/main.cpp
)

# Download Box2D
function(download_box2d)
    custom_download_project(Box2D-cmake
        GIT_REPOSITORY https://github.com/tobanteGaming/Box2D-cmake.git
        GIT_TAG        6cf2fe06fbd7a5b6f37affcbe39a4f3a7063df85
    )
endfunction()

# Box2D
if(NOT TARGET tobanteGaming::Box2D)
    download_box2d()
    add_subdirectory(${FIXING_COLLISIONS_EXTERNAL}/Box2D-cmake box2d)
endif()
target_link_libraries(Box2D-comparison PUBLIC tobanteGaming::Box2D)

# Eigen through libigl
target_link_libraries(Box2D-comparison PUBLIC igl::core)

# Warnings
target_link_libraries(Box2D-comparison PUBLIC warnings::all)

# Json Parser
target_link_libraries(Box2D-comparison PUBLIC nlohmann_json::nlohmann_json)

# Logger
target_link_libraries(Box2D-comparison PUBLIC spdlog::spdlog)

# Command-Line Interface
target_link_libraries(Box2D-comparison PUBLIC CLI11::CLI11)
