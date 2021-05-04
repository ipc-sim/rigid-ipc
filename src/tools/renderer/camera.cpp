#include "camera.hpp"

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <igl/PI.h>
#include <igl/barycenter.h>

#include "eigen_json.hpp"

namespace swr {

Camera::Camera(const nlohmann::json& json)
    : zoom(1)
    , shift(0, 0, 0)
{
    is_perspective = json["is_perspective"];
    from_json(json["position"], position);
    from_json(json["gaze"], gaze);
    from_json(json["view_up"], view_up);
    field_of_view = json["field_of_view"];
    field_of_view =
        std::min(Float(172.847), std::max(Float(0.367), field_of_view));
    field_of_view *= igl::PI / 180;
    orthographic_scale = json["orthographic_scale"];
    near_clip = json["near_clip"];
    far_clip = json["far_clip"];
    from_json(json["resolution"], resolution);
}

Matrix4F Camera::shift_and_zoom() const
{
    return (Eigen::Scaling(zoom) * Eigen::Translation<Float, 3>(shift))
        .matrix();
}

Matrix4F Camera::view_transform() const
{
    Vector3F w = -gaze.normalized();
    Vector3F u = view_up.cross(w).normalized();
    Vector3F v = w.cross(u);

    Matrix4F view_inv;
    // clang-format off
    view_inv << u.x(), v.x(), w.x(), position.x(),
                u.y(), v.y(), w.y(), position.y(),
                u.z(), v.z(), w.z(), position.z(),
                0, 0, 0, 1;
    // clang-format on

    return view_inv.inverse();
}

Matrix4F Camera::projection_matrix() const
{
    Float t;
    if (is_perspective) {
        t = tan(field_of_view / 2) * near_clip;
    } else {
        t = orthographic_scale / 2;
    }
    Float b = -t;
    Float aspect_ratio = float(resolution.x()) / resolution.y();
    Float r = t * aspect_ratio;
    Float l = -r;
    Float n = -near_clip;
    Float f = -far_clip;

    Matrix4F M_orth;
    // clang-format off
    M_orth << 2 / (r - l), 0, 0, -(r + l) / (r - l),
              0, 2 / (t - b), 0, -(t + b) / (t - b),
              0, 0, 2 / (n - f), -(n + f) / (n - f),
              0, 0, 0, 1;
    // clang-format on

    Matrix4F P;
    if (is_perspective) {
        // clang-format off
        P << n, 0, 0, 0,
             0, n, 0, 0,
             0, 0, n + f, -f * n,
             0, 0, 1, 0;
        // clang-format on
    } else {
        P.setIdentity();
    }

    return M_orth * P;
}

void get_scale_and_shift_to_fit_mesh(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    Float& zoom,
    Vector3F& shift)
{
    if (V.rows() == 0) {
        return;
    }

    Eigen::MatrixXd BC;
    if (F.rows() <= 1) {
        BC = V;
    } else {
        igl::barycenter(V, F, BC);
    }

    auto min_point = V.colwise().minCoeff();
    auto max_point = V.colwise().maxCoeff();
    auto centroid = (0.5 * (min_point + max_point)).eval();
    shift.setConstant(0);
    shift.head(centroid.size()) = -centroid.cast<Float>();
    zoom = 2.0 / (max_point - min_point).array().abs().maxCoeff();
}

void Camera::align_camera_center(
    const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    if (V.rows() == 0) {
        return;
    }
    get_scale_and_shift_to_fit_mesh(V, F, zoom, shift);
}

} // namespace swr
