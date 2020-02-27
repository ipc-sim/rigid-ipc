#include "distance_barrier_constraint.hpp"

#include <igl/slice_mask.h>

#include <ccd/rigid_body_collision_detection.hpp>
#include <ccd/time_of_impact.hpp>
#include <geometry/distance.hpp>

#include <logger.hpp>
#include <profiler.hpp>

namespace ccd {
namespace opt {

    NLOHMANN_JSON_SERIALIZE_ENUM(
        BarrierType,
        { { BarrierType::IPC, "ipc" },
          { BarrierType::POLY_LOG, "poly_log" },
          { BarrierType::SPLINE, "spline" } })

    DistanceBarrierConstraint::DistanceBarrierConstraint()
        : DistanceBarrierConstraint("distance_barrier_constraint")
    {
    }

    DistanceBarrierConstraint::DistanceBarrierConstraint(
        const std::string& name)
        : CollisionConstraint(name)
        , custom_inital_epsilon(1.0)
        , min_distance(1E-10)
        , active_constraint_scale(1.5)
        , barrier_type(BarrierType::POLY_LOG)
        , m_barrier_epsilon(0.0)
    {
    }

    void DistanceBarrierConstraint::settings(const nlohmann::json& json)
    {
        CollisionConstraint::settings(json);
        custom_inital_epsilon = json["custom_initial_epsilon"].get<double>();
        active_constraint_scale = json["active_constraint_scale"].get<double>();
        min_distance = json["min_distance"].get<double>();
        barrier_type = json["barrier_type"].get<BarrierType>();
    }

    nlohmann::json DistanceBarrierConstraint::settings() const
    {
        nlohmann::json json = CollisionConstraint::settings();
        json["custom_inital_epsilon"] = custom_inital_epsilon;
        json["active_constraint_scale"] = active_constraint_scale;
        json["min_distance"] = min_distance;
        json["barrier_type"] = barrier_type;

        return json;
    }

    void DistanceBarrierConstraint::initialize()
    {
        m_barrier_epsilon = custom_inital_epsilon;
        CollisionConstraint::initialize();
    }

    bool DistanceBarrierConstraint::has_active_collisions(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1) const
    {
        Candidates candidates;

        detect_collision_candidates(
            bodies, poses_t0, poses_t1 - poses_t0,
            dim_to_collision_type(bodies.dim()), candidates, detection_method,
            /*inflation_radius=*/0);

        for (const auto& ev_candidate : candidates.ev_candidates) {
            double toi, alpha;
            bool are_colliding = detect_edge_vertex_collisions_narrow_phase(
                bodies, poses_t0, poses_t1 - poses_t0, ev_candidate, toi,
                alpha);
            if (are_colliding) {
                return true;
            }
        }
        for (const auto& ee_candidate : candidates.ee_candidates) {
            double toi, edge0_alpha, edge1_alpha;
            bool are_colliding = detect_edge_edge_collisions_narrow_phase(
                bodies, poses_t0, poses_t1 - poses_t0, ee_candidate, toi,
                edge0_alpha, edge1_alpha);
            if (are_colliding) {
                return true;
            }
        }
        for (const auto& fv_candidate : candidates.fv_candidates) {
            double toi, u, v;
            bool are_colliding = detect_face_vertex_collisions_narrow_phase(
                bodies, poses_t0, poses_t1 - poses_t0, fv_candidate, toi, u, v);
            if (are_colliding) {
                return true;
            }
        }
        return false;
    }

    void DistanceBarrierConstraint::compute_constraints(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses,
        const physics::Poses<double>& displacements,
        Eigen::VectorXd& barriers)
    {
        Candidates candidates;
        construct_active_barrier_set(bodies, poses, displacements, candidates);
        compute_candidates_constraints(
            bodies, poses, displacements, candidates, barriers);
    }

    void DistanceBarrierConstraint::construct_active_barrier_set(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses,
        const physics::Poses<double>& displacements,
        Candidates& barriers) const
    {
        NAMED_PROFILE_POINT(
            "distance_barrier__update_collision_set", BROAD_PHASE)
        NAMED_PROFILE_POINT("distance_barrier__update_active_set", NARROW_PHASE)

        PROFILE_START(BROAD_PHASE)

        Candidates candidates;

        detect_collision_candidates(
            bodies, poses, displacements, dim_to_collision_type(bodies.dim()),
            candidates, detection_method,
            /*inflation_radius=*/active_constraint_scale * m_barrier_epsilon);

        PROFILE_END(BROAD_PHASE)

        PROFILE_START(NARROW_PHASE)

        barriers.clear();

        // Compute the world vertices at the end of the time-step
        Eigen::MatrixXd vertices_t1 =
            bodies.world_vertices(poses + displacements);

        // Compute the edge-vertex distance candidates
        for (const auto& ev_candidate : candidates.ev_candidates) {
            double distance = ccd::geometry::point_segment_distance<double>(
                vertices_t1.row(ev_candidate.vertex_index),
                vertices_t1.row(bodies.m_edges(ev_candidate.edge_index, 0)),
                vertices_t1.row(bodies.m_edges(ev_candidate.edge_index, 1)));

            if (distance < active_constraint_scale * m_barrier_epsilon) {
                barriers.ev_candidates.push_back(ev_candidate);
            }
        }

        // TODO: 3D
        if (bodies.dim() != 2) {
            throw NotImplementedError(
                "construct_active_barrier_set() not implmented for 3D!");
        }

        PROFILE_END(NARROW_PHASE)
    }

    void DistanceBarrierConstraint::debug_compute_distances(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses,
        const physics::Poses<double>& displacements,
        Eigen::VectorXd& distances) const
    {
        Candidates candidates;
        construct_active_barrier_set(bodies, poses, displacements, candidates);

        typedef double T;

        Eigen::MatrixX<T> vertices_t1 =
            bodies.world_vertices(poses + displacements);

        distances.resize(candidates.size(), 1);
        distances.setConstant(T(0.0));
        for (size_t i = 0; i < candidates.ev_candidates.size(); i++) {
            const auto& ev_candidate = candidates.ev_candidates[i];
            // a and b are the endpoints of the edge; c is the vertex
            long edge_id = ev_candidate.edge_index;
            int a_id = bodies.m_edges(edge_id, 0);
            int b_id = bodies.m_edges(edge_id, 1);
            long c_id = ev_candidate.vertex_index;
            // Check that the vertex is not an endpoint of the edge
            assert(a_id != c_id && b_id != c_id);
            Eigen::VectorX3<T> a = vertices_t1.row(a_id);
            Eigen::VectorX3<T> b = vertices_t1.row(b_id);
            Eigen::VectorX3<T> c = vertices_t1.row(c_id);

            distances(int(i)) =
                ccd::geometry::point_segment_distance<T>(c, a, b);
        }

        if (bodies.dim() != 2) {
            throw NotImplementedError(
                "DistanceBarrierConstraint::debug_compute_distances() not "
                "implemented in 3D!");
        }
    }

    double DistanceBarrierConstraint::distance_barrier_grad(
        const double distance, const double eps)
    {
        return barrier_gradient(distance, eps, barrier_type);
    }

    Eigen::VectorXd DistanceBarrierConstraint::distance_barrier_grad(
        const Eigen::VectorXd& a,
        const Eigen::VectorXd& b,
        const Eigen::VectorXd& c)
    {
        // TODO: 3D
        // 6 dof = 3 points * 2 dof / point
        typedef AutodiffType<6> Diff;
        Diff::activate();
        Diff::D1VectorXd da = Diff::d1vars(0, a);
        Diff::D1VectorXd db = Diff::d1vars(2, b);
        Diff::D1VectorXd dc = Diff::d1vars(4, c);
        Diff::DDouble1 barrier = distance_barrier<Diff::DDouble1>(da, db, dc);

        return barrier.getGradient();
    }

    Eigen::MatrixXd DistanceBarrierConstraint::distance_barrier_hess(
        const Eigen::VectorXd& a,
        const Eigen::VectorXd& b,
        const Eigen::VectorXd& c)
    {
        // TODO: 3D
        // 6 dof = 3 points * 2 dof / point
        typedef AutodiffType<6> Diff;
        Diff::activate();
        Diff::D2VectorXd da = Diff::d2vars(0, a);
        Diff::D2VectorXd db = Diff::d2vars(2, b);
        Diff::D2VectorXd dc = Diff::d2vars(4, c);
        Diff::DDouble2 barrier = distance_barrier<Diff::DDouble2>(da, db, dc);
        return barrier.getHessian();
    }

} // namespace opt
} // namespace ccd
