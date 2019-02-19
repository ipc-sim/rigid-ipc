#include <algorithm> // std::sort
#include <iostream>

#include "state.hpp"

#include <ccd/collision_volume.hpp>
#include <ccd/prune_impacts.hpp>

#include <io/read_scene.hpp>
#include <io/write_scene.hpp>

namespace ccd {
State::State()
    : canvas_width(10)
    , canvas_height(10)
    , current_impact(-1)
{
}

void State::load_scene(std::string filename)
{
    io::read_scene(filename, vertices, edges, displacements);
    volumes.resize(edges.rows());
    volumes.setZero();
    volume_grad.resize(vertices.rows(), kDIM);
    volume_grad.setZero();

    current_impact = -1;
    time = 0.0;
    selected_displacements.clear();
    selected_points.clear();
    impacts = nullptr;
}

void State::save_scene(std::string filename)
{
    io::write_scene(filename, vertices, edges, displacements);
}
// CRUD Scene ---------------------------------------------------
void State::add_vertex(const Eigen::RowVector2d& position)
{
    long lastid = vertices.rows();
    vertices.conservativeResize(lastid + 1, kDIM);
    vertices.row(lastid) << position;

    displacements.conservativeResize(lastid + 1, kDIM);
    displacements.row(lastid) << 0.0, -0.1;

    volume_grad.conservativeResize(lastid + 1, kDIM);
    volume_grad.row(lastid).setConstant(0.0);
}

void State::add_edges(const Eigen::MatrixX2i& new_edges)
{
    assert(new_edges.cols() == 2);

    long lastid = edges.rows();
    edges.conservativeResize(lastid + new_edges.rows(), 2);
    for (unsigned i = 0; i < new_edges.rows(); ++i)
        edges.row(lastid + i) << new_edges.row(i);

    // Add a new rows to the volume vector
    volumes.conservativeResize(lastid + new_edges.rows());
    // Add zero volumes for each new edge
    for (unsigned i = 0; i < new_edges.rows(); ++i)
        volumes(lastid + i) = 0.0;
}

void State::set_vertex_position(
    const int vertex_idx, const Eigen::RowVector2d& position)
{
    vertices.row(vertex_idx) = position;
}

void State::move_vertex(const int vertex_idx, const Eigen::RowVector2d& delta)
{
    vertices.row(vertex_idx) += delta;
}
void State::move_displacement(
    const int vertex_idx, const Eigen::RowVector2d& delta)
{
    displacements.row(vertex_idx) += delta;
}

// CCD ---------------------------------------------------------
void State::detect_edge_vertex_collisions()
{
    impacts = ccd::detect_edge_vertex_collisions(
        vertices, displacements, edges, detection_method);
    std::sort(impacts->begin(), impacts->end(),
        EdgeVertexImpact::compare_impacts_by_time);
}

void State::compute_collision_volumes()
{
    EdgeEdgeImpactsPtr ee_impacts;
    ee_impacts = ccd::EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
        this->edges, this->impacts);

    auto pruned_impacts = ccd::prune_impacts(ee_impacts);
    std::cout << "# of Impacts: " << pruned_impacts->size() << std::endl;

    for (auto impact : *pruned_impacts) {
        int edge_id = impact.first;
        EdgeEdgeImpactPtr ee_impact = impact.second;

        this->volumes(edge_id) = ccd::collision_volume(
            vertices, displacements, edges, ee_impact, this->epsilon);

        // TODO: Add gradient of volume computation
    }

    std::cout << "Collision Volumes:\n" << (this->volumes) << std::endl;
}


void State::run_full_pipeline()
{
    this->detect_edge_vertex_collisions();
    this->compute_collision_volumes();
}

}
