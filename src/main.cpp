#include <igl/opengl/glfw/Viewer.h>

#include <logger.hpp>
#include <viewer/UISimState.hpp>

int main(int argc, char* argv[])
{
    spdlog::set_level(spdlog::level::info);
    ccd::UISimState ui;
    ui.launch();
}
