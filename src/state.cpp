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
    , current_ev_impact(-1)
    , current_ee_impact(-1)
    , min_edge_width(0.0)
{
}

void State::load_scene(std::string filename)
{
    io::read_scene(filename, vertices, edges, displacements);

    reset_scene();
}

void State::reset_scene()
{
    volumes.resize(edges.rows());
    volumes.setZero();

    edges_impact.resize(edges.rows());
    edges_impact.setConstant(nullptr);

    volume_grad.resize(0, 0);
    volume_grad.setZero();

    current_ee_impact = -1;
    current_ev_impact = -1;
    time = 0.0;
    selected_displacements.clear();
    selected_points.clear();
    ev_impacts = nullptr;
}

void State::save_scene(std::string filename)
{
    io::write_scene(filename, vertices, edges, displacements);
}
// -----------------------------------------------------------------------------
// CRUD Scene
// -----------------------------------------------------------------------------
void State::add_vertex(const Eigen::RowVector2d& position)
{
    long lastid = vertices.rows();
    vertices.conservativeResize(lastid + 1, kDIM);
    vertices.row(lastid) << position;

    displacements.conservativeResize(lastid + 1, kDIM);
    displacements.row(lastid) << 0.0, -0.1;

    reset_impacts();
}

void State::add_edges(const Eigen::MatrixX2i& new_edges)
{
    assert(new_edges.cols() == 2);

    long lastid = edges.rows();
    edges.conservativeResize(lastid + new_edges.rows(), 2);
    for (unsigned i = 0; i < new_edges.rows(); ++i)
        edges.row(lastid + i) << new_edges.row(i);

    // Add a new rows to the volume vector
    edges_impact.conservativeResize(lastid + new_edges.rows());
    volumes.conservativeResize(lastid + new_edges.rows());

    reset_impacts();
}

void State::set_vertex_position(
    const int vertex_idx, const Eigen::RowVector2d& position)
{
    vertices.row(vertex_idx) = position;
    reset_impacts();
}

void State::move_vertex(const int vertex_idx, const Eigen::RowVector2d& delta)
{
    vertices.row(vertex_idx) += delta;
    reset_impacts();
}

void State::move_displacement(
    const int vertex_idx, const Eigen::RowVector2d& delta)
{
    displacements.row(vertex_idx) += delta;
    reset_impacts();
}

Eigen::MatrixX2d State::get_vertex_at_time(){
    return vertices + displacements * double(time);;
}

Eigen::MatrixX2d State::get_volume_grad(){
    if (current_ee_impact < 0 || volume_grad.cols() == 0){
        return Eigen::MatrixX2d::Zero(vertices.rows(), kDIM);
    }
    Eigen::MatrixXd grad = volume_grad.col(current_ee_impact);
    grad.resize( grad.rows() / kDIM, kDIM);

    return grad;
}


// -----------------------------------------------------------------------------
// CCD
// -----------------------------------------------------------------------------
void State::reset_impacts()
{

    volumes.setZero();
    volume_grad.resize(0, 0);
    edges_impact.setConstant(nullptr);

    if (ev_impacts != nullptr) {
        ev_impacts->clear();
    }
    if (ee_impacts != nullptr) {
        ee_impacts->clear();
    }
}

void State::detect_edge_vertex_collisions()
{
    // get impacts between vertex and edge
    ev_impacts = ccd::detect_edge_vertex_collisions(
        vertices, displacements, edges, detection_method);

    std::sort(ev_impacts->begin(), ev_impacts->end(),
        EdgeVertexImpact::compare_impacts_by_time);

    // transform to impacts between two edges
    EdgeEdgeImpactsPtr ee_impacts;
    ee_impacts = ccd::EdgeEdgeImpact::convert_edge_vertex_to_edge_edge_impacts(
        this->edges, this->ev_impacts);

    // assign first impact to each edge
    auto pruned_impacts = ccd::prune_impacts(ee_impacts);
    std::cout << "# of EE-Impacts: " << pruned_impacts->size() << std::endl;

    // we store one impact for each edge on edges_impact
    // and the impacts in ee_impacts;
    this->ee_impacts = std::make_shared<EdgeEdgeImpacts>();
    for (auto impact : *pruned_impacts) {
        int edge_id = impact.first;
        EdgeEdgeImpactPtr ee_impact = impact.second;
        this->edges_impact[edge_id] = ee_impact;
        this->ee_impacts->push_back(ee_impact);
    }

    std::sort(ee_impacts->begin(), ee_impacts->end(),
        EdgeEdgeImpact::compare_impacts_by_time);
}

void State::compute_collision_volumes()
{

    size_t num_impacts = this->ee_impacts->size();
    this->volume_grad.resize(kDIM * this->vertices.rows(), long(num_impacts));

    for (long i = 0; i < this->edges_impact.rows(); ++i) {
        EdgeEdgeImpactPtr ee_impact = this->edges_impact[i];

        if (ee_impact == nullptr) {
            this->volumes(i) = 0;
            continue;
        }
        // get collision volume for this edge
        this->volumes(i) = ccd::collision_volume(
            vertices, displacements, edges, ee_impact, int(i), this->epsilon);

        // get collision volume gradient for this edge
        Eigen::VectorXd grad;
        ccd::collision_volume_grad(
            vertices, displacements, edges, ee_impact, int(i), this->epsilon, grad);
        assert(grad.rows() == this->volume_grad.rows());
        this->volume_grad.col(i) = grad;
    }

    std::cout << "Collision Volumes:\n"
              << (this->volumes) << std::endl;
}

void State::run_full_pipeline()
{
    this->detect_edge_vertex_collisions();
    this->compute_collision_volumes();
}

}
