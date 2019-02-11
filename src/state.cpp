#include <algorithm> // std::sort

#include "state.hpp"

#include "autogen/auto_test_function.hpp"
#include <FixingCollisions/collision_volume.hpp>
#include <read_scene.hpp>
#include <write_scene.hpp>
<<<<<<< HEAD
=======
#include <FixingCollisions/collision_volume.hpp>
#include "autogen/test_function.hpp"
>>>>>>> 1861a0e2fcb3751bf76848fd29f6048075b49e94

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
<<<<<<< HEAD
void State::detect_edge_vertex_collisions()
{
    impacts = ccd::detect_edge_vertex_collisions(
        vertices, displacements, edges, detection_method);
    std::sort(impacts->begin(), impacts->end(),
        EdgeVertexImpact::compare_impacts_by_time);
    // NOTE: dummy calls just so the code is called
    for (size_t i = 0; i < impacts->size(); i++) {
        ccd::collision_volume(vertices, displacements, edges, *impacts->at(i),
            epsilon, Eigen::VectorXd());
        auto edge = edges.row(impacts->at(i)->edge_index);
        ccd::autogen::test_function(vertices.row(edge(0)),
            vertices.row(edge(1)), vertices.row(edge(0)), vertices.row(edge(1)),
            displacements.row(edge(0)), displacements.row(edge(1)),
            displacements.row(edge(0)), displacements.row(edge(1)), epsilon);
    }
=======
void State::detect_edge_vertex_collisions(){
    impacts = ccd::detect_edge_vertex_collisions(vertices, displacements, edges, detection_method);
>>>>>>> 1861a0e2fcb3751bf76848fd29f6048075b49e94
}

}
