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

    DistanceBarrierConstraint::DistanceBarrierConstraint(
        const std::string& name)
        : CollisionConstraint(name)
        , initial_barrier_activation_distance(1e-2)
        , min_distance(0.0)
        , active_constraint_scale(1.01)
        , barrier_type(BarrierType::POLY_LOG)
        , m_barrier_activation_distance(0.0)
    {
    }

    void DistanceBarrierConstraint::settings(const nlohmann::json& json)
    {
        CollisionConstraint::settings(json);
        initial_barrier_activation_distance =
            json["initial_barrier_activation_distance"].get<double>();
        active_constraint_scale = json["active_constraint_scale"].get<double>();
        min_distance = json["min_distance"].get<double>();
        barrier_type = json["barrier_type"].get<BarrierType>();
    }

    nlohmann::json DistanceBarrierConstraint::settings() const
    {
        nlohmann::json json = CollisionConstraint::settings();
        json["initial_barrier_activation_distance"] =
            initial_barrier_activation_distance;
        json["active_constraint_scale"] = active_constraint_scale;
        json["min_distance"] = min_distance;
        json["barrier_type"] = barrier_type;
        return json;
    }

    void DistanceBarrierConstraint::initialize()
    {
        m_barrier_activation_distance = initial_barrier_activation_distance;
        CollisionConstraint::initialize();
    }

    bool DistanceBarrierConstraint::has_active_collisions(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1) const
    {
        PROFILE_POINT("collisions_detection");
        NAMED_PROFILE_POINT("collisions_detection__narrow_phase", NARROW_PHASE);

        // This function will profile itself
        Candidates candidates;
        detect_collision_candidates(
            bodies, poses_t0, poses_t1, dim_to_collision_type(bodies.dim()),
            candidates, detection_method, trajectory_type,
            /*inflation_radius=*/0);

        PROFILE_START();
        PROFILE_START(NARROW_PHASE)
        bool has_collisions = has_active_collisions_narrow_phase(
            bodies, poses_t0, poses_t1, candidates);
        PROFILE_END(NARROW_PHASE)
        PROFILE_END();

        return has_collisions;
    }

    bool DistanceBarrierConstraint::has_active_collisions_narrow_phase(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1,
        const Candidates& candidates) const
    {
        for (const auto& ev_candidate : candidates.ev_candidates) {
            double toi, alpha;
            bool are_colliding = detect_edge_vertex_collisions_narrow_phase(
                bodies, poses_t0, poses_t1, ev_candidate, toi, alpha,
                trajectory_type);
            if (are_colliding) {
                return true;
            }
        }
        for (const auto& ee_candidate : candidates.ee_candidates) {
            double toi, edge0_alpha, edge1_alpha;
            bool are_colliding = detect_edge_edge_collisions_narrow_phase(
                bodies, poses_t0, poses_t1, ee_candidate, toi, edge0_alpha,
                edge1_alpha, trajectory_type);
            if (are_colliding) {
                return true;
            }
        }
        for (const auto& fv_candidate : candidates.fv_candidates) {
            double toi, u, v;
            bool are_colliding = detect_face_vertex_collisions_narrow_phase(
                bodies, poses_t0, poses_t1, fv_candidate, toi, u, v,
                trajectory_type);
            if (are_colliding) {
                return true;
            }
        }
        return false;
    }

    void DistanceBarrierConstraint::construct_active_barrier_set(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses,
        Candidates& barriers) const
    {
        NAMED_PROFILE_POINT(
            "distance_barrier__update_collision_set", BROAD_PHASE)
        NAMED_PROFILE_POINT("distance_barrier__update_active_set", NARROW_PHASE)

        PROFILE_START(BROAD_PHASE)

        Candidates candidates;

        // Use discrete collision detection
        // There is no trajectory, so use the more efficient linearized version
        // TODO: Replace all of these duplicate V or poses with a single
        //       parameter version
        detect_collision_candidates(
            bodies, poses, poses, dim_to_collision_type(bodies.dim()),
            candidates, detection_method, TrajectoryType::LINEARIZED,
            /*inflation_radius=*/active_constraint_scale
                * m_barrier_activation_distance);

        PROFILE_END(BROAD_PHASE)

        PROFILE_START(NARROW_PHASE)

        barriers.clear();

        // Compute the distance candidates
        // TODO: Consider skipping this step
        // TODO: Store these distance for use later
        double activation_distance =
            active_constraint_scale * m_barrier_activation_distance;
        tbb::parallel_invoke(
            [&] {
                // TODO: Consider parallelizing this loop
                for (const auto& ev_candidate : candidates.ev_candidates) {
                    double distance =
                        ccd::geometry::point_segment_distance<double>(
                            bodies.world_vertex(
                                poses, ev_candidate.vertex_index),
                            bodies.world_vertex(
                                poses,
                                bodies.m_edges(ev_candidate.edge_index, 0)),
                            bodies.world_vertex(
                                poses,
                                bodies.m_edges(ev_candidate.edge_index, 1)));

                    if (distance < activation_distance) {
                        barriers.ev_candidates.push_back(ev_candidate);
                    }
                }
            },
            [&] {
                for (const auto& ee_candidate : candidates.ee_candidates) {
                    double distance =
                        ccd::geometry::segment_segment_distance<double>(
                            bodies.world_vertex(
                                poses,
                                bodies.m_edges(ee_candidate.edge0_index, 0)),
                            bodies.world_vertex(
                                poses,
                                bodies.m_edges(ee_candidate.edge0_index, 1)),
                            bodies.world_vertex(
                                poses,
                                bodies.m_edges(ee_candidate.edge1_index, 0)),
                            bodies.world_vertex(
                                poses,
                                bodies.m_edges(ee_candidate.edge1_index, 1)));

                    if (distance < activation_distance) {
                        barriers.ee_candidates.push_back(ee_candidate);
                    }
                }
            },
            [&] {
                for (const auto& fv_candidate : candidates.fv_candidates) {
                    double distance =
                        ccd::geometry::point_triangle_distance<double>(
                            bodies.world_vertex(
                                poses, fv_candidate.vertex_index),
                            bodies.world_vertex(
                                poses,
                                bodies.m_faces(fv_candidate.face_index, 0)),
                            bodies.world_vertex(
                                poses,
                                bodies.m_faces(fv_candidate.face_index, 1)),
                            bodies.world_vertex(
                                poses,
                                bodies.m_faces(fv_candidate.face_index, 2)));

                    if (distance < activation_distance) {
                        barriers.fv_candidates.push_back(fv_candidate);
                    }
                }
            });

        PROFILE_END(NARROW_PHASE)
    }

    void DistanceBarrierConstraint::compute_distances(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses,
        Eigen::VectorXd& distances) const
    {
        Candidates candidates;
        construct_active_barrier_set(bodies, poses, candidates);

        distances.resize(candidates.size());

        size_t offset = 0;
        for (size_t i = 0; i < candidates.ev_candidates.size(); i++) {
            const auto& ev_candidate = candidates.ev_candidates[i];
            distances(i + offset) =
                ccd::geometry::point_segment_distance<double>(
                    bodies.world_vertex(poses, ev_candidate.vertex_index),
                    bodies.world_vertex(
                        poses, bodies.m_edges(ev_candidate.edge_index, 0)),
                    bodies.world_vertex(
                        poses, bodies.m_edges(ev_candidate.edge_index, 1)));
        }

        offset += candidates.ev_candidates.size();
        for (size_t i = 0; i < candidates.ee_candidates.size(); i++) {
            const auto& ee_candidate = candidates.ee_candidates[i];
            distances(i + offset) =
                ccd::geometry::segment_segment_distance<double>(
                    bodies.world_vertex(
                        poses, bodies.m_edges(ee_candidate.edge0_index, 0)),
                    bodies.world_vertex(
                        poses, bodies.m_edges(ee_candidate.edge0_index, 1)),
                    bodies.world_vertex(
                        poses, bodies.m_edges(ee_candidate.edge1_index, 0)),
                    bodies.world_vertex(
                        poses, bodies.m_edges(ee_candidate.edge1_index, 1)));
        }

        offset += candidates.ee_candidates.size();
        for (size_t i = 0; i < candidates.fv_candidates.size(); i++) {
            const auto& fv_candidate = candidates.fv_candidates[i];
            distances(i + offset) =
                ccd::geometry::point_triangle_distance<double>(
                    bodies.world_vertex(poses, fv_candidate.vertex_index),
                    bodies.world_vertex(
                        poses, bodies.m_faces(fv_candidate.face_index, 0)),
                    bodies.world_vertex(
                        poses, bodies.m_faces(fv_candidate.face_index, 1)),
                    bodies.world_vertex(
                        poses, bodies.m_faces(fv_candidate.face_index, 2)));
        }
    }

} // namespace opt
} // namespace ccd
