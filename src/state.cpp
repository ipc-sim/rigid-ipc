#include "state.hpp"

namespace ccd {
State::State()
    : canvas_width(10)
    , canvas_height(10)
{
}

void State::load_scene(std::string /*filename*/)
{
    int num_vertices = 3;
    int num_edges = 1;

    // surface
    vertices.resize(num_vertices, 3);
    vertices.row(0) << -1.0, 0.0, 0.0;
    vertices.row(1) << 1.0, 0.0, 0.0;
    vertices.row(2) << 0.0, 0.5, 0.0;

    edges.resize(num_edges, 2);
    edges.row(0) << 0, 1;

    // per-node displacements
    displacements.resize(num_vertices, 3);
    displacements.row(0) << 0.0, 0.5, 0.0;
    displacements.row(1) << 0.0, 0.5, 0.0;
    displacements.row(2) << 0.0, -1.0, 0.0;

    // per-node volume gradient ?
    volume_grad.resize(num_vertices, 3);
    volume_grad.setZero();
}

void State::add_vertex(const Eigen::RowVector3d& position){
    long lastid = vertices.rows();
    vertices.conservativeResize(lastid + 1, 3);
    vertices.row(lastid) << position;

    displacements.conservativeResize(lastid + 1, 3);
    displacements.row(lastid) << 0.0, -0.1, 0.0;

    volume_grad.conservativeResize(lastid + 1, 3);
    volume_grad.row(lastid).setConstant(0.0);

}

void State::add_edges(const Eigen::MatrixXi& new_edges)
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

}
