#pragma once

#include <nlohmann/json.hpp>

#include "float.hpp"

namespace swr {

class Camera {
public:
    Camera() = default;
    Camera(const nlohmann::json& json);

    Matrix4F shift_and_zoom() const;
    Matrix4F view_transform() const;
    Matrix4F projection_matrix() const;

    void
    align_camera_center(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

    Vector3F position;
    Vector3F gaze;
    Vector3F view_up;

    bool is_perspective;
    Float field_of_view; // between 0 and PI
    Float orthographic_scale;
    Float near_clip;
    Float far_clip;
    Eigen::Vector2i resolution;

    Float zoom;
    Vector3F shift;
};

} // namespace swr
