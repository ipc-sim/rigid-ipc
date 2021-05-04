#include "render_mesh.hpp"

#include <iostream>

#include <Eigen/Geometry>
#include <Eigen/LU>

#include "eigen_json.hpp"
#include "raster.hpp"
#include "write_png.hpp"

namespace swr {

Scene::Scene(const nlohmann::json& scene_json)
{
    // Read scene info
    from_json(scene_json["background_color"], background_color);
    from_json(scene_json["ambient_light"], ambient_light);
    line_thickness = scene_json["line_thickness"];

    // Read camera info
    camera = Camera(scene_json["camera"]);

    // Read materials
    materials.emplace_back(scene_json["static_body_material"]);    // #B3B3B3
    materials.emplace_back(scene_json["kinematic_body_material"]); // #FF8000
    materials.emplace_back(scene_json["dynamic_body_material"]);   // #E74C3C

    // Read lights
    for (const auto& entry : scene_json["lights"]) {
        lights.emplace_back(entry);
    }
}

UniformAttributes build_uniform(const Scene& scene)
{
    UniformAttributes uniform;
    uniform.lights = scene.lights;
    uniform.camera = scene.camera;
    uniform.ambient = scene.ambient_light;
    const Camera& camera = uniform.camera;
    uniform.ST = camera.shift_and_zoom();
    uniform.M =
        camera.projection_matrix() * camera.view_transform() * uniform.ST;
    return uniform;
}

Shaders create_face_shaders()
{
    Shaders shaders;

    shaders.vertex_shader = [](const VertexAttributes& va,
                               const UniformAttributes& uniform) {
        VertexAttributes out;
        out.position = uniform.M * va.position;
        out.position /= out.position.w();
        Vector3F ambient_color =
            va.material.ambient_color.cwiseProduct(uniform.ambient);

        // Skip this step if the diffuse and specular are not used
        Vector3F lights_color(0, 0, 0);
        if (va.material.diffuse_color.squaredNorm() != 0
            || va.material.specular_color.squaredNorm() != 0) {

            // The vertex position has to be after translation and scaling
            Vector4F homogenous_hit = uniform.ST * va.position;
            Vector3F hit = homogenous_hit.head<3>() / homogenous_hit.w();
            Vector3F N = va.normal;
            Vector3F V;
            if (uniform.camera.is_perspective) {
                V = (uniform.camera.position - hit).normalized();
            } else {
                V = -uniform.camera.gaze.normalized();
            }
            for (auto light : uniform.lights) {
                Vector3F Li = (light.position - hit).normalized();

                // Diffuse contribution
                Vector3F diffuse =
                    va.material.diffuse_color * std::max(Li.dot(N), Float(0));

                // Specular contribution
                Vector3F H = (Li + V).normalized();
                Vector3F specular = va.material.specular_color
                    * std::pow(std::max(N.dot(H), Float(0)),
                               va.material.specular_exponent);

                // Attenuate lights according to the squared distance to the
                // lights
                Vector3F D = light.position - hit;
                lights_color +=
                    (diffuse + specular).cwiseProduct(light.intensity)
                    / D.squaredNorm();
            }
        }

        out.color = ambient_color + lights_color;
        return out;
    };

    shaders.fragment_shader = [](const VertexAttributes& va,
                                 const UniformAttributes& uniform) {
        FragmentAttributes out(va.color);
        // out.color.head<3>() = va.bc;
        // if (va.bc.x() < 0.01 || va.bc.y() < 0.01 || va.bc.z() < 0.01) {
        //     out.color << 0.0, 0.0, 0.0, 1.0;
        // }
        out.depth = -va.position(2);
        return out;
    };

    shaders.blending_shader = [](const FragmentAttributes& fa,
                                 FrameBufferAttributes& pixel) {
        if (fa.depth < pixel.depth) {
            pixel.color =
                round(fa.color.array().max(0).min(1) * 255).cast<uint8_t>();
            pixel.depth = fa.depth;
        }
    };

    return shaders;
}

Shaders create_edge_shaders()
{
    Shaders shaders;

    shaders.vertex_shader = [](const VertexAttributes& va,
                               const UniformAttributes& uniform) {
        VertexAttributes out;
        out.position = uniform.M * va.position;
        out.position /= out.position.w();
        return out;
    };

    shaders.fragment_shader = [](const VertexAttributes& va,
                                 const UniformAttributes& uniform) {
        FragmentAttributes out(0, 0, 0);
        out.depth = -va.position(2);
        return out;
    };

    shaders.blending_shader = [](const FragmentAttributes& fa,
                                 FrameBufferAttributes& pixel) {
        if (fa.depth < pixel.depth + 1e-4) // offset
        {
            pixel.color =
                round(fa.color.array().max(0).min(1) * 255).cast<uint8_t>();
            // out.depth = fa.depth;
        }
    };

    return shaders;
}

// Render a mesh with vertices V, edges E, faces F, and vertex colors C
bool render_mesh(
    const Scene& scene,
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& E,
    const Eigen::MatrixXi& F,
    const Eigen::VectorXi& C,
    const std::string& filename)
{
    UniformAttributes uniform = build_uniform(scene);

    FrameBuffer frame_buffer = create_frame_buffer(
        scene.camera.resolution,
        round(scene.background_color.array().max(0).min(1) * 255)
            .cast<uint8_t>());

    std::vector<VertexAttributes> face_vertex_attributes;
    for (int i = 0; i < F.rows(); i++) {
        Vector3F l1 = (V.row(F(i, 1)) - V.row(F(i, 0))).cast<Float>();
        Vector3F l2 = (V.row(F(i, 2)) - V.row(F(i, 0))).cast<Float>();
        Vector3F normal = (l1).cross(l2).normalized();
        for (int j = 0; j < 3; j++) {
            VertexAttributes va(V.row(F(i, j)).cast<Float>());
            va.material = scene.materials[C[F(i, j)]];
            va.normal = normal;
            face_vertex_attributes.push_back(va);
        }
    }
    rasterize_triangles(
        create_face_shaders(), uniform, face_vertex_attributes, frame_buffer);

    if (scene.line_thickness > 0) {
        // draw wireframe on top
        std::vector<VertexAttributes> edge_vertex_attributes;
        for (int i = 0; i < E.rows(); i++) {
            for (int j = 0; j < E.cols(); j++) {
                VertexAttributes va(V.row(E(i, j)).cast<Float>());
                va.material = scene.materials[C[E(i, j)]];
                va.normal.setZero();
                edge_vertex_attributes.push_back(va);
            }
        }
        rasterize_lines(
            create_edge_shaders(), uniform, edge_vertex_attributes,
            scene.line_thickness, frame_buffer);
    }

    // Convert Framebuffer to RGBA matrices
    Eigen::Matrix<uint8_t, Eigen::Dynamic, Eigen::Dynamic> R, G, B, A;
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

    return write_png(R, G, B, A, filename);
}

} // namespace swr
