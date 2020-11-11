#pragma once

#include <Eigen/Core>
#include <string>
#include <vector>

#include "attributes.hpp"

namespace swr {

// Contains the three shaders used by the rasterizer
class Shaders {
public:
    // Vertex Shader
    std::function<VertexAttributes(
        const VertexAttributes&, const UniformAttributes&)>
        vertex_shader;
    // Fragment Shader
    std::function<FragmentAttributes(
        const VertexAttributes&, const UniformAttributes&)>
        fragment_shader;
    // Blending Shader
    std::function<void(const FragmentAttributes&, FrameBufferAttributes&)>
        blending_shader;
};

// Rasterizes a single triangle v1,v2,v3 using the provided shaders and
// uniforms. Note: v1, v2, and v3 needs to be in the canonical view volume (i.e.
// after being processed by the vertex shader)
void rasterize_triangle(
    const Shaders& shaders,
    const UniformAttributes& uniform,
    const VertexAttributes& v1,
    const VertexAttributes& v2,
    const VertexAttributes& v3,
    FrameBuffer& frameBuffer);

// Rasterizes a collection of triangles, assembling one triangle for each 3
// consecutive vertices. Note: the vertices will be processed by the vertex
// shader
void rasterize_triangles(
    const Shaders& shaders,
    const UniformAttributes& uniform,
    const std::vector<VertexAttributes>& vertices,
    FrameBuffer& frameBuffer);

// Rasterizes a single line v1,v2 of thickness line_thickness using the provided
// shaders and uniforms. Note: v1, v2 needs to be in the canonical view volume
// (i.e. after being processed by the vertex shader)
void rasterize_line(
    const Shaders& shaders,
    const UniformAttributes& uniform,
    const VertexAttributes& v1,
    const VertexAttributes& v2,
    Float line_thickness,
    FrameBuffer& frameBuffer);

// Rasterizes a collection of lines, assembling one line for each 2 consecutive
// vertices. Note: the vertices will be processed by the vertex shader
void rasterize_lines(
    const Shaders& shaders,
    const UniformAttributes& uniform,
    const std::vector<VertexAttributes>& vertices,
    Float line_thickness,
    FrameBuffer& frameBuffer);

// Exports the framebuffer to a uint8 raw image
void framebuffer_to_uint8(
    const FrameBuffer& frameBuffer, std::vector<uint8_t>& image);

} // namespace swr
