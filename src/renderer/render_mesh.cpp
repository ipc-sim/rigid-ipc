#include "render_mesh.hpp"

#include <iostream>

#include <Eigen/Geometry>
#include <Eigen/LU>
#include <fmt/format.h>
#include <igl/PI.h>
#include <igl/png/writePNG.h>
#include <spdlog/spdlog.h>

#include <renderer/raster.hpp>

Scene::Scene(const nlohmann::json& scene_json)
{
    // Helper function to read a Vector3F from a json array
    auto read_vec2i = [](const nlohmann::json& x) {
        return Eigen::Vector2i(x[0], x[1]);
    };
    auto read_vec3 = [](const nlohmann::json& x) {
        return Eigen::Vector3F(x[0], x[1], x[2]);
    };
    auto read_vec4 = [](const nlohmann::json& x) {
        return Eigen::Vector4F(x[0], x[1], x[2], x[3]);
    };

    // Read scene info
    background_color = read_vec4(scene_json["background_color"]);
    ambient_light = read_vec3(scene_json["ambient_light"]);
    line_thickness = scene_json["line_thickness"];

    // Read camera info
    camera.is_perspective = scene_json["camera"]["is_perspective"];
    camera.position = read_vec3(scene_json["camera"]["position"]);
    camera.gaze = read_vec3(scene_json["camera"]["gaze"]);
    camera.view_up = read_vec3(scene_json["camera"]["view_up"]);
    camera.field_of_view = scene_json["camera"]["field_of_view"];
    camera.field_of_view =
        std::min(Float(172.847), std::max(Float(0.367), camera.field_of_view));
    camera.field_of_view *= igl::PI / 180;
    camera.orthographic_scale = scene_json["camera"]["orthographic_scale"];
    camera.near_clip = scene_json["camera"]["near_clip"];
    camera.far_clip = scene_json["camera"]["far_clip"];
    camera.resolution = read_vec2i(scene_json["camera"]["resolution"]);

    // Read materials
    for (const auto& entry : scene_json["materials"]) {
        Material mat;
        mat.ambient_color = read_vec3(entry["ambient"]);
        mat.diffuse_color = read_vec3(entry["diffuse"]);
        mat.specular_color = read_vec3(entry["specular"]);
        mat.specular_exponent = entry["shininess"];
        materials.push_back(mat);
    }

    // Read lights
    for (const auto& entry : scene_json["lights"]) {
        Light light;
        light.position = read_vec3(entry["position"]);
        light.intensity = read_vec3(entry["color"]);
        lights.push_back(light);
    }
}

UniformAttributes build_uniform(const Scene& scene)
{
    const Camera& camera = scene.camera;

    UniformAttributes uniform;

    uniform.lights = scene.lights;
    uniform.camera = scene.camera;
    uniform.ambient = scene.ambient_light;

    Eigen::Vector3F w = -camera.gaze.normalized();
    Eigen::Vector3F u = camera.view_up.cross(w).normalized();
    Eigen::Vector3F v = w.cross(u);

    Eigen::Matrix4F M_cam_inv;
    // clang-format off
	M_cam_inv << u(0), v(0), w(0), camera.position(0),
                 u(1), v(1), w(1), camera.position(1),
                 u(2), v(2), w(2), camera.position(2),
                 0, 0, 0, 1;
    // clang-format on
    Eigen::Matrix4F M_cam = M_cam_inv.inverse();

    Float t;
    if (camera.is_perspective) {
        t = tan(camera.field_of_view / 2) * camera.near_clip;
    } else {
        t = camera.orthographic_scale / 2;
    }
    Float b = -t;
    Float aspect_ratio = float(camera.resolution.x()) / camera.resolution.y();
    Float r = t * aspect_ratio;
    Float l = -r;
    Float n = -camera.near_clip;
    Float f = -camera.far_clip;

    Eigen::Matrix4F M_orth;
    // clang-format off
    M_orth << 2 / (r - l), 0, 0, -(r + l) / (r - l),
              0, 2 / (t - b), 0, -(t + b) / (t - b),
              0, 0, 2 / (n - f), -(n + f) / (n - f),
              0, 0, 0, 1;
    // clang-format on

    Eigen::Matrix4F P;
    if (camera.is_perspective) {
        // clang-format off
        P << n, 0, 0, 0,
             0, n, 0, 0,
             0, 0, n + f, -f * n,
             0, 0, 1, 0;
        // clang-format on
    } else {
        P.setIdentity();
    }

    uniform.M = M_orth * P * M_cam;

    return uniform;
}

Program create_shading_program()
{
    Program program;

    program.vertex_shader = [](const VertexAttributes& va,
                               const UniformAttributes& uniform) {
        VertexAttributes out;
        out.position = uniform.M * va.position;
        out.position /= out.position.w();
        Eigen::Vector3F color =
            va.material.ambient_color.cwiseProduct(uniform.ambient);

        // Eigen::Vector3F hit(va.position(0), va.position(1), va.position(2));
        // for (auto light : uniform.lights) {
        //     Eigen::Vector3F Li = (light.position - hit).normalized();
        //     Eigen::Vector3F N = va.normal;
        //     Eigen::Vector3F diffuse =
        //         va.material.diffuse_color * std::max(Li.dot(N), Float(0));
        //     Eigen::Vector3F H;
        //     if (uniform.camera.is_perspective) {
        //         H = (Li + (uniform.camera.position - hit).normalized())
        //                 .normalized();
        //     } else {
        //         H = (Li - uniform.camera.gaze.normalized()).normalized();
        //     }
        //     Eigen::Vector3F specular = va.material.specular_color
        //         * std::pow(std::max(N.dot(H), Float(0)),
        //                    va.material.specular_exponent);
        //     Eigen::Vector3F D = light.position - hit;
        //     color += (diffuse + specular).cwiseProduct(light.intensity)
        //         / D.squaredNorm();
        // }
        out.color = color;
        return out;
    };

    program.fragment_shader = [](const VertexAttributes& va,
                                 const UniformAttributes& uniform) {
        FragmentAttributes out(va.color(0), va.color(1), va.color(2));
        out.depth = -va.position(2);
        return out;
    };

    program.blending_shader = [](const FragmentAttributes& fa,
                                 FrameBufferAttributes& pixel) {
        if (fa.depth < pixel.depth) {
            pixel.color = (fa.color * 255).cast<uint8_t>();
            pixel.depth = fa.depth;
        }
    };

    return program;
}

// Render a mesh with vertices V, edges E, faces F, and vertex colors C
bool render_mesh(
    const Scene& scene,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Eigen::MatrixXd& C,
    const std::string& filename)
{
    Program program = create_shading_program();

    UniformAttributes uniform = build_uniform(scene);

    FrameBuffer frame_buffer(
        scene.camera.resolution.x(), scene.camera.resolution.y());
    for (int i = 0; i < frame_buffer.rows(); i++) {
        for (int j = 0; j < frame_buffer.cols(); j++) {
            frame_buffer(i, j).color =
                (scene.background_color * 255).cast<uint8_t>();
        }
    }

    std::vector<VertexAttributes> face_vertex_attributes;
    for (int i = 0; i < F.rows(); i++) {
        Eigen::Vector3F l1 = (V.row(F(i, 1)) - V.row(F(i, 0))).cast<Float>();
        Eigen::Vector3F l2 = (V.row(F(i, 2)) - V.row(F(i, 0))).cast<Float>();
        Eigen::Vector3F normal = (l1).cross(l2).normalized();
        for (int j = 0; j < 3; j++) {
            VertexAttributes va(V.row(F(i, j)).cast<Float>());
            va.material.ambient_color = C.row(F(i, j)).cast<Float>();
            va.material.diffuse_color << 0.5, 0.5, 0.5;
            va.material.specular_color.setZero();
            va.material.specular_exponent = 0;
            va.normal = normal;
            face_vertex_attributes.push_back(va);
        }
    }
    rasterize_triangles(program, uniform, face_vertex_attributes, frame_buffer);

    // shaders for wireframe
    Program program_wire;
    program_wire.vertex_shader = [](const VertexAttributes& va,
                                    const UniformAttributes& uniform) {
        VertexAttributes out;
        out.position = uniform.M * va.position;
        out.position /= out.position.w();
        return out;
    };

    program_wire.fragment_shader = [](const VertexAttributes& va,
                                      const UniformAttributes& uniform) {
        FragmentAttributes out(0, 0, 0);
        out.depth = -va.position(2);
        return out;
    };

    program_wire.blending_shader = [](const FragmentAttributes& fa,
                                      FrameBufferAttributes& pixel) {
        if (fa.depth < pixel.depth + 1e-4) // offset
        {
            pixel.color = (fa.color * 255).cast<uint8_t>();
            // out.depth = fa.depth;
        }
    };

    // draw wireframe on top
    std::vector<VertexAttributes> edge_vertex_attributes;
    for (int i = 0; i < E.rows(); i++) {
        for (int j = 0; j < E.cols(); j++) {
            VertexAttributes va(V.row(E(i, j)).cast<Float>());
            va.material.ambient_color.setZero();
            va.material.diffuse_color.setZero();
            va.material.specular_color.setZero();
            va.material.specular_exponent = 0;
            va.normal.setZero();
            edge_vertex_attributes.push_back(va);
        }
    }
    rasterize_lines(
        program_wire, uniform, edge_vertex_attributes, scene.line_thickness,
        frame_buffer);

    // Convert Framebuffer to RGBA matrices
    Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> R, G, B, A;
    R.resizeLike(frame_buffer);
    G.resizeLike(frame_buffer);
    B.resizeLike(frame_buffer);
    A.resizeLike(frame_buffer);
    for (int i = 0; i < frame_buffer.rows(); i++) {
        for (int j = 0; j < frame_buffer.cols(); j++) {
            R(i, j) = frame_buffer(i, j).color(0);
            G(i, j) = frame_buffer(i, j).color(1);
            B(i, j) = frame_buffer(i, j).color(2);
            A(i, j) = frame_buffer(i, j).color(3);
        }
    }

    return igl::png::writePNG(R, G, B, A, filename);
}
