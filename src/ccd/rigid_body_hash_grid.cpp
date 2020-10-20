// A spatial hash grid for rigid bodies with angular trajectories.
#include "rigid_body_hash_grid.hpp"

#include <cmath>
#include <iostream>

#include <interval/interval.hpp>
#include <logger.hpp>

namespace ccd {

void compute_vertices_intervals(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    Eigen::MatrixX<Interval>& vertices)
{
    Interval t(0, 1);

    typedef physics::Poses<Interval> PosesI;

    PosesI posesI_t0 = physics::cast<double, Interval>(poses_t0);
    PosesI posesI_t1 = physics::cast<double, Interval>(poses_t1);

    PosesI posesI = physics::interpolate(posesI_t0, posesI_t1, t);

    vertices = bodies.world_vertices(posesI);
}

void RigidBodyHashGrid::resize(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const double inflation_radius)
{
    Eigen::MatrixX<Interval> vertices;
    compute_vertices_intervals(bodies, poses_t0, poses_t1, vertices);

    // Compute the extents of the mesh
    Eigen::VectorX3<Interval> mesh_extents(bodies.dim());
    for (int i = 0; i < vertices.rows(); i++) {
        for (int j = 0; j < vertices.cols(); j++) {
            mesh_extents(j) = i == 0
                ? vertices(i, j)
                : boost::numeric::hull(mesh_extents(j), vertices(i, j));
        }
    }
    Eigen::VectorX3d min(mesh_extents.size()), max(mesh_extents.size());
    for (int i = 0; i < mesh_extents.size(); i++) {
        min(i) = mesh_extents(i).lower() - inflation_radius;
        max(i) = mesh_extents(i).upper() + inflation_radius;
    }

    // Compute the average displacement length
    double average_displacement_length = 0;
    for (int i = 0; i < vertices.rows(); i++) {
        double max_side_width = 0;
        for (int j = 0; j < vertices.cols(); j++) {
            max_side_width =
                std::max(max_side_width, boost::numeric::width(vertices(i, j)));
        }
        average_displacement_length += max_side_width;
    }

    average_displacement_length /= vertices.rows();
    double average_edge_length = bodies.average_edge_length;

    double cell_size =
        std::max(average_displacement_length, average_edge_length)
        + inflation_radius;
    HashGrid::resize(min, max, cell_size);
}

void RigidBodyHashGrid::addBodies(
    const physics::RigidBodyAssembler& bodies,
    const physics::Poses<double>& poses_t0,
    const physics::Poses<double>& poses_t1,
    const double inflation_radius)
{
    assert(bodies.num_bodies() == poses_t0.size());
    Eigen::MatrixX<Interval> vertices;
    compute_vertices_intervals(bodies, poses_t0, poses_t1, vertices);

    // Create a bounding box for all vertices
    tbb::concurrent_vector<AABB> vertices_aabb;
    Eigen::VectorX3d min(bodies.dim()), max(bodies.dim());
    vertices_aabb.reserve(vertices.rows());
    for (long i = 0; i < vertices.rows(); i++) {
        for (int j = 0; j < vertices.cols(); j++) {
            min(j) = vertices(i, j).lower() - inflation_radius;
            max(j) = vertices(i, j).upper() + inflation_radius;
        }
        vertices_aabb.emplace_back(min, max);
    }

    // Add vertice, edges, and faces in parallel
    tbb::parallel_invoke(
        [&] {
            // Add all vertices of the bodies
            for (long i = 0; i < vertices.rows(); i++) {
                this->addElement(vertices_aabb[i], i, this->m_vertexItems);
            }
        },
        [&] {
            // Add all edge of the bodies
            for (long i = 0; i < bodies.m_edges.rows(); i++) {
                this->addElement(
                    AABB(
                        vertices_aabb[bodies.m_edges(i, 0)],
                        vertices_aabb[bodies.m_edges(i, 1)]),
                    i, this->m_edgeItems);
            }
        },
        [&] {
            // Add all faces of the bodies
            for (long i = 0; i < bodies.m_faces.rows(); i++) {
                this->addElement(
                    AABB(
                        AABB(
                            vertices_aabb[bodies.m_faces(i, 0)],
                            vertices_aabb[bodies.m_faces(i, 1)]),
                        vertices_aabb[bodies.m_faces(i, 2)]),
                    i, this->m_faceItems);
            }
        });
}

} // namespace ccd
