#pragma once

#include <array>
#include <vector>

#include <Eigen/Core>
#include <tbb/mutex.h>

using Float = double;
namespace Eigen {
typedef Matrix<Float, 3, 3> Matrix3F;
typedef Matrix<Float, 4, 4> Matrix4F;
typedef Matrix<Float, 2, 1> Vector2F;
typedef Matrix<Float, 3, 1> Vector3F;
typedef Matrix<Float, 4, 1> Vector4F;
} // namespace Eigen

struct Light {
    Eigen::Vector3F position;
    Eigen::Vector3F intensity;
};

struct Material {
    Eigen::Vector3F ambient_color;
    Eigen::Vector3F diffuse_color;
    Eigen::Vector3F specular_color;
    Float specular_exponent; // Also called "shininess"
};

struct Camera {
    bool is_perspective;
    Eigen::Vector3F position;
    Eigen::Vector3F gaze;
    Eigen::Vector3F view_up;
    Float field_of_view; // between 0 and PI
    Float orthographic_scale;
    Float near_clip;
    Float far_clip;
    Eigen::Vector2i resolution;
};

class VertexAttributes {
public:
    VertexAttributes(Float x = 0, Float y = 0, Float z = 0, Float w = 1)
        : VertexAttributes(Eigen::Vector3F(x, y, z), w)
    {
    }
    VertexAttributes(const Eigen::Vector3F& pos, Float w = 1)
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
        return r;
    }

    Eigen::Vector4F position;
    Eigen::Vector3F normal;
    Eigen::Vector3F color;
    Material material;
};

class FragmentAttributes {
public:
    FragmentAttributes(Float r = 0, Float g = 0, Float b = 0, Float a = 1)
    {
        color << r, g, b, a;
    }
    Float depth;
    Eigen::Vector4F color;
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
    tbb::mutex lock;
};

class UniformAttributes {
public:
    Eigen::Matrix4F M;
    std::vector<Light> lights;
    Eigen::Vector3F ambient;
    Camera camera;
};
