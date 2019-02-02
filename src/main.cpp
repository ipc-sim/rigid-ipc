#include <igl/opengl/glfw/Viewer.h>

#include "viewer.hpp"

int main(int /*argc*/, char* /*argv*/[])
{

    igl::opengl::glfw::Viewer viewer;
    ccd::ViewerMenu menu;
    viewer.plugins.push_back(&menu);

    viewer.core.set_rotation_type(igl::opengl::ViewerCore::ROTATION_TYPE_NO_ROTATION);
    viewer.core.orthographic = true;
    viewer.core.is_animating = true;
    viewer.core.lighting_factor = 0.0;

    viewer.launch();
}
