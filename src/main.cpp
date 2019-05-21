#include <igl/opengl/glfw/Viewer.h>

#include <viewer/viewer.hpp>
#include <logger.hpp>

int main(int argc, char* argv[])
{
    spdlog::set_level(spdlog::level::info);

    igl::opengl::glfw::Viewer viewer;
    ccd::ViewerMenu menu;

    std::string scene_file;
    if (argc > 1) {
        menu.scene_file = argv[1];
    }

    viewer.plugins.push_back(&menu);

    viewer.core().set_rotation_type(
        igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION);
    viewer.core().orthographic = true;
    viewer.core().is_animating = true;
    viewer.core().lighting_factor = 0.0;

    viewer.launch();
}
