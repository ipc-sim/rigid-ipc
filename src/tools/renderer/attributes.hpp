#pragma once

#include <array>
#include <iostream>
#include <mutex>
#include <vector>

#include <Eigen/Core>

#include "camera.hpp"
#include "eigen_json.hpp"
#include "float.hpp"

namespace swr {

struct Light {
    Light() = default;
    Light(const nlohmann::json& json)
    {
        from_json(json["position"], position);
        from_json(json["color"], intensity);
    }

    Vector3F position;
    Vector3F intensity;
};

struct Material {
    Material() = default;
    Material(const nlohmann::json& json)
    {
        from_json(json["ambient"], ambient_color);
        from_json(json["diffuse"], diffuse_color);
        from_json(json["specular"], specular_color);
        specular_exponent = json["shininess"];
    }

    Vector3F ambient_color;
    Vector3F diffuse_color;
    Vector3F specular_color;
    Float specular_exponent; // Also called "shininess"
};

class VertexAttributes {
public:
    VertexAttributes(Float x = 0, Float y = 0, Float z = 0, Float w = 1)
        : VertexAttributes(Vector3F(x, y, z), w)
    {
    }
    VertexAttributes(const Vector3F& pos, Float w = 1)
    {
        position.head<3>() = pos;
        position.w() = w;
    }

    // Interpolates the vertex attributes
    static VertexAttributes interpolate(
        const VertexAttributes& a,
        const VertexAttributes& b,
        const VertexAttributes& c,
        const Float alpha,
        const Float beta,
        const Float gamma)
    {
        VertexAttributes r;
        r.position =
            alpha * a.position + beta * b.position + gamma * c.position;
        r.color = alpha * a.color + beta * b.color + gamma * c.color;
        r.bc << alpha, beta, gamma;
        return r;
    }

    Vector4F position;
    Vector3F normal;
    Vector3F color;
    Vector3F bc;
    Material material;
};

class FragmentAttributes {
public:
    FragmentAttributes(Float r = 0, Float g = 0, Float b = 0, Float a = 1)
        : FragmentAttributes(Vector3F(r, g, b), a)
    {
    }
    FragmentAttributes(const Vector3F& rgb, Float a = 1)
    {
        color.head<3>() = rgb;
        color(3) = a;
    }
    Float depth;
    Vector4F color;
};

class FrameBufferAttributes {
public:
    FrameBufferAttributes(
        uint8_t r = 0, uint8_t g = 0, uint8_t b = 0, uint8_t a = 255)
        : FrameBufferAttributes(Eigen::Matrix<uint8_t, 4, 1>(r, g, b, a))
    {
    }
    FrameBufferAttributes(const Eigen::Matrix<uint8_t, 4, 1>& color)
        : color(color)
        , depth(std::numeric_limits<Float>::infinity())
    {
    }
    Float depth;
    Eigen::Matrix<uint8_t, 4, 1> color;
    std::mutex lock;
};

typedef Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic>
    FrameBuffer;

inline FrameBuffer create_frame_buffer(
    const Eigen::Vector2i& resolution,
    const Eigen::Matrix<uint8_t, 4, 1>& bg_color)
{
    FrameBuffer frame_buffer(resolution.x(), resolution.y());
    for (int i = 0; i < frame_buffer.rows(); i++) {
        for (int j = 0; j < frame_buffer.cols(); j++) {
            frame_buffer(i, j).color = bg_color;
        }
    }
    return frame_buffer;
}

class UniformAttributes {
public:
    Matrix4F M;
    Matrix4F ST;
    std::vector<Light> lights;
    Vector3F ambient;
    Camera camera;
};

} // namespace swr
