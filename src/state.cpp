#include "state.hpp"

#include <read_scene.hpp>
#include <write_scene.hpp>

namespace ccd {
State::State()
    : canvas_width(10)
    , canvas_height(10),
      current_impact(-1)
{
}

void State::load_scene(std::string filename)
{
    io::read_scene(filename, vertices, edges, displacements);
    volume_grad.resize(vertices.rows(), 3);
    volume_grad.setZero();
}

void State::save_scene(std::string filename)
{
    io::write_scene(filename, vertices, edges, displacements);
}
// CRUD Scene ---------------------------------------------------
void State::add_vertex(const Eigen::RowVector3d& position)
{
    long lastid = vertices.rows();
    vertices.conservativeResize(lastid + 1, 3);
    vertices.row(lastid) << position;

    displacements.conservativeResize(lastid + 1, 3);
    displacements.row(lastid) << 0.0, -0.1, 0.0;

    volume_grad.conservativeResize(lastid + 1, 3);
    volume_grad.row(lastid).setConstant(0.0);
}

void State::add_edges(const Eigen::MatrixX2i& new_edges)
{
    assert(new_edges.cols() == 2);

    long lastid = edges.rows();
    edges.conservativeResize(lastid + new_edges.rows(), 2);
    for (unsigned i = 0; i < new_edges.rows(); ++i)
        edges.row(lastid + i) << new_edges.row(i);
}

void State::set_vertex_position(const int vertex_idx, const Eigen::RowVector3d& position)
{
    vertices.row(vertex_idx) = position;
}

void State::move_vertex(const int vertex_idx, const Eigen::RowVector3d& delta)
{
    vertices.row(vertex_idx) += delta;
}
void State::move_displacement(const int vertex_idx, const Eigen::RowVector3d& delta)
{
    displacements.row(vertex_idx) += delta;
}

// CCD ---------------------------------------------------------
void State::detect_edge_vertex_collisions(){
    impacts = ccd::detect_edge_vertex_collisions(vertices, vertices + displacements, edges, detection_method);
}

}
