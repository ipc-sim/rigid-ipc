#include <igl/opengl/glfw/Viewer.h>

#include <logger.hpp>
#include <viewer/UISimState.hpp>

int main(int argc, char* argv[])
{
    // tbb::task_scheduler_init init(1);
    spdlog::set_level(spdlog::level::info);
    ccd::UISimState ui;
    ui.launch();
}
