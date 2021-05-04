#include "raster.hpp"

#include <iostream>

#include <Eigen/LU> // Needed for .inverse()

#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>

namespace swr {

void rasterize_triangle(
    const Shaders& shaders,
    const UniformAttributes& uniform,
    const VertexAttributes& v1,
    const VertexAttributes& v2,
    const VertexAttributes& v3,
    FrameBuffer& frameBuffer)
{
    // Collect coordinates into a matrix and convert to canonical representation
    Eigen::Matrix<Float, 3, 4> p;
    p.row(0) = v1.position.array() / v1.position.w();
    p.row(1) = v2.position.array() / v2.position.w();
    p.row(2) = v3.position.array() / v3.position.w();

    // Coordinates are in -1..1, rescale to pixel size (x,y only)
    p.col(0) = ((p.col(0).array() + 1.0) / 2.0) * frameBuffer.rows();
    p.col(1) = ((p.col(1).array() + 1.0) / 2.0) * frameBuffer.cols();

    // Find bounding box in pixels
    int lx = std::floor(p.col(0).minCoeff());
    int ly = std::floor(p.col(1).minCoeff());
    int ux = std::ceil(p.col(0).maxCoeff());
    int uy = std::ceil(p.col(1).maxCoeff());

    // Clamp to framebuffer
    lx = std::min(std::max(lx, int(0)), int(frameBuffer.rows() - 1));
    ly = std::min(std::max(ly, int(0)), int(frameBuffer.cols() - 1));
    ux = std::max(std::min(ux, int(frameBuffer.rows() - 1)), int(0));
    uy = std::max(std::min(uy, int(frameBuffer.cols() - 1)), int(0));

    // Build the implicit triangle representation
    Matrix3F A;
    A.col(0) = p.row(0).head<3>();
    A.col(1) = p.row(1).head<3>();
    A.col(2) = p.row(2).head<3>();
    A.row(2) << 1.0, 1.0, 1.0;

    Matrix3F Ai = A.inverse();

    // Rasterize the triangle
    tbb::parallel_for(
        tbb::blocked_range2d<size_t>(lx, ux + 1, ly, uy + 1),
        [&](const tbb::blocked_range2d<size_t>& r) {
            for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
                for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
                    // The pixel center is offset by 0.5, 0.5
                    Vector3F pixel(i + 0.5, j + 0.5, 1);
                    Vector3F b = Ai * pixel;
                    if (b.minCoeff() >= 0) {
                        VertexAttributes va = VertexAttributes::interpolate(
                            v1, v2, v3, b[0], b[1], b[2]);
                        // Only render fragments within the bi-unit cube
                        if (va.position.z() >= -1 && va.position.z() <= 1) {
                            FragmentAttributes frag =
                                shaders.fragment_shader(va, uniform);
                            frameBuffer(i, j).lock.lock();
                            shaders.blending_shader(frag, frameBuffer(i, j));
                            frameBuffer(i, j).lock.unlock();
                        }
                    }
                }
            }
        });
}

void rasterize_triangles(
    const Shaders& shaders,
    const UniformAttributes& uniform,
    const std::vector<VertexAttributes>& vertices,
    FrameBuffer& frameBuffer)
{
    // Call vertex shader on all vertices (parallel)
    std::vector<VertexAttributes> v(vertices.size());
    tbb::parallel_for(size_t(0), vertices.size(), [&](size_t i) {
        v[i] = shaders.vertex_shader(vertices[i], uniform);
    });

    // Call the rasterization function on every triangle
    assert(vertices.size() % 3 == 0);
    tbb::parallel_for(size_t(0), vertices.size() / 3, [&](size_t i) {
        rasterize_triangle(
            shaders, uniform, v[i * 3 + 0], v[i * 3 + 1], v[i * 3 + 2],
            frameBuffer);
    });
}

void rasterize_line(
    const Shaders& shaders,
    const UniformAttributes& uniform,
    const VertexAttributes& v1,
    const VertexAttributes& v2,
    Float line_thickness,
    FrameBuffer& frameBuffer)
{
    // Collect coordinates into a matrix and convert to canonical
    // representation
    Eigen::Matrix<Float, 2, 4> p;
    p.row(0) = v1.position.array() / v1.position[3];
    p.row(1) = v2.position.array() / v2.position[3];

    // Coordinates are in -1..1, rescale to pixel size (x,y only)
    p.col(0) = ((p.col(0).array() + 1.0) / 2.0) * frameBuffer.rows();
    p.col(1) = ((p.col(1).array() + 1.0) / 2.0) * frameBuffer.cols();

    // Find bounding box in pixels, adding the line thickness
    int lx = std::floor(p.col(0).minCoeff() - line_thickness);
    int ly = std::floor(p.col(1).minCoeff() - line_thickness);
    int ux = std::ceil(p.col(0).maxCoeff() + line_thickness);
    int uy = std::ceil(p.col(1).maxCoeff() + line_thickness);

    // Clamp to framebuffer
    lx = std::min(std::max(lx, int(0)), int(frameBuffer.rows() - 1));
    ly = std::min(std::max(ly, int(0)), int(frameBuffer.cols() - 1));
    ux = std::max(std::min(ux, int(frameBuffer.rows() - 1)), int(0));
    uy = std::max(std::min(uy, int(frameBuffer.cols() - 1)), int(0));

    // We only need the 2d coordinates of the endpoints of the line
    Vector2F l1(p(0, 0), p(0, 1));
    Vector2F l2(p(1, 0), p(1, 1));

    // Parametrize the line as l1 + t (l2-l1)
    Float t = -1;
    Float ll = (l1 - l2).squaredNorm();

    // Rasterize the line
    tbb::parallel_for(
        tbb::blocked_range2d<size_t>(lx, ux + 1, ly, uy + 1),
        [&](const tbb::blocked_range2d<size_t>& r) {
            for (size_t i = r.rows().begin(); i != r.rows().end(); ++i) {
                for (size_t j = r.cols().begin(); j != r.cols().end(); ++j) {
                    // The pixel center is offset by 0.5, 0.5
                    Vector2F pixel(i + 0.5, j + 0.5);

                    if (ll == 0.0)
                        // The segment has zero length
                        t = 0;
                    else {
                        // Project p on the line
                        t = (pixel - l1).dot(l2 - l1) / ll;
                        // Clamp between 0 and 1
                        t = std::fmax(0, std::fmin(1, t));
                    }

                    Vector2F pixel_p = l1 + t * (l2 - l1);

                    if ((pixel - pixel_p).squaredNorm()
                        < (line_thickness * line_thickness)) {
                        VertexAttributes va = VertexAttributes::interpolate(
                            v1, v2, v1, 1 - t, t, 0);
                        // Only render fragments within the bi-unit cube
                        if (va.position[2] >= -1 && va.position[2] <= 1) {
                            FragmentAttributes frag =
                                shaders.fragment_shader(va, uniform);
                            frameBuffer(i, j).lock.lock();
                            shaders.blending_shader(frag, frameBuffer(i, j));
                            frameBuffer(i, j).lock.unlock();
                        }
                    }
                }
            }
        });
}

void rasterize_lines(
    const Shaders& shaders,
    const UniformAttributes& uniform,
    const std::vector<VertexAttributes>& vertices,
    Float line_thickness,
    FrameBuffer& frameBuffer)
{
    // Call vertex shader on all vertices (parallel)
    std::vector<VertexAttributes> v(vertices.size());
    tbb::parallel_for(size_t(0), vertices.size(), [&](size_t i) {
        v[i] = shaders.vertex_shader(vertices[i], uniform);
    });

    // Call the rasterization function on every line
    assert(vertices.size() % 2 == 0);
    tbb::parallel_for(size_t(0), vertices.size() / 2, [&](size_t i) {
        rasterize_line(
            shaders, uniform, v[i * 2 + 0], v[i * 2 + 1], line_thickness,
            frameBuffer);
    });
}

void framebuffer_to_uint8(
    const FrameBuffer& frameBuffer, std::vector<uint8_t>& image)
{
    const int w = frameBuffer.rows();     // Image width
    const int h = frameBuffer.cols();     // Image height
    const int comp = 4;                   // 4 Channels Red, Green, Blue, Alpha
    const int stride_in_bytes = w * comp; // Length of one row in bytes
    image.resize(w * h * comp, 0);        // The image itself;

    for (unsigned wi = 0; wi < w; ++wi) {
        for (unsigned hi = 0; hi < h; ++hi) {
            unsigned hif = h - 1 - hi;
            image[(hi * w * 4) + (wi * 4) + 0] = frameBuffer(wi, hif).color[0];
            image[(hi * w * 4) + (wi * 4) + 1] = frameBuffer(wi, hif).color[1];
            image[(hi * w * 4) + (wi * 4) + 2] = frameBuffer(wi, hif).color[2];
            image[(hi * w * 4) + (wi * 4) + 3] = frameBuffer(wi, hif).color[3];
        }
    }
}

} // namespace swr
