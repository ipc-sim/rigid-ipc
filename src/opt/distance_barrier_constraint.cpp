#include "distance_barrier_constraint.hpp"

#include <mutex>
#include <tbb/parallel_for_each.h>

#include <igl/slice_mask.h>
#include <ipc/ipc.hpp>

#include <ccd/rigid/broad_phase.hpp>
#include <ccd/rigid/rigid_body_hash_grid.hpp>
#include <geometry/distance.hpp>
#include <io/serialize_json.hpp>
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
        , barrier_type(BarrierType::IPC)
        , m_barrier_activation_distance(0.0)
    {
    }

    void DistanceBarrierConstraint::settings(const nlohmann::json& json)
    {
        CollisionConstraint::settings(json);
        initial_barrier_activation_distance =
            json["initial_barrier_activation_distance"].get<double>();
        barrier_type = json["barrier_type"].get<BarrierType>();
    }

    nlohmann::json DistanceBarrierConstraint::settings() const
    {
        nlohmann::json json = CollisionConstraint::settings();
        json["initial_barrier_activation_distance"] =
            initial_barrier_activation_distance;
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
        PROFILE_POINT("DistanceBarrierConstraint::has_active_collisions");
        NAMED_PROFILE_POINT(
            "DistanceBarrierConstraint::has_active_collisions_narrow_phase",
            NARROW_PHASE);

        PROFILE_START();
        // This function will profile itself
        ipc::Candidates candidates;
        detect_collision_candidates(
            bodies, poses_t0, poses_t1, dim_to_collision_type(bodies.dim()),
            candidates, detection_method, trajectory_type,
            /*inflation_radius=*/0);

        PROFILE_START(NARROW_PHASE)
        bool has_collisions = has_active_collisions_narrow_phase(
            bodies, poses_t0, poses_t1, candidates);
        PROFILE_END(NARROW_PHASE)
        PROFILE_END();

        return has_collisions;
    }

    static int counter = 0;

    void save_candidate(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1,
        const ipc::EdgeVertexCandidate& ev_candidate)
    {
        nlohmann::json json;
        long bodyA_id, vertex_id, bodyB_id, edge_id;
        bodies.global_to_local_vertex(
            ev_candidate.vertex_index, bodyA_id, vertex_id);
        bodies.global_to_local_edge(ev_candidate.edge_index, bodyB_id, edge_id);
        const auto& bodyA = bodies[bodyA_id];
        const auto& bodyB = bodies[bodyB_id];
        json["type"] = "ev";

        json["edge"] = nlohmann::json();
        json["edge"]["vertex0"] =
            io::to_json(bodyB.vertices.row(bodyB.edges(edge_id, 0)));
        json["edge"]["vertex1"] =
            io::to_json(bodyB.vertices.row(bodyB.edges(edge_id, 1)));
        json["edge"]["pose_t0"] = nlohmann::json();
        json["edge"]["pose_t0"]["position"] =
            io::to_json(poses_t0[bodyB_id].position);
        json["edge"]["pose_t0"]["rotation"] =
            io::to_json(poses_t0[bodyB_id].rotation);
        json["edge"]["pose_t1"] = nlohmann::json();
        json["edge"]["pose_t1"]["position"] =
            io::to_json(poses_t1[bodyB_id].position);
        json["edge"]["pose_t1"]["rotation"] =
            io::to_json(poses_t1[bodyB_id].rotation);

        json["vertex"] = nlohmann::json();
        json["vertex"]["vertex"] = io::to_json(bodyA.vertices.row(vertex_id));
        json["vertex"]["pose_t0"]["position"] =
            io::to_json(poses_t0[bodyA_id].position);
        json["vertex"]["pose_t0"]["rotation"] =
            io::to_json(poses_t0[bodyA_id].rotation);
        json["vertex"]["pose_t1"] = nlohmann::json();
        json["vertex"]["pose_t1"]["position"] =
            io::to_json(poses_t1[bodyA_id].position);
        json["vertex"]["pose_t1"]["rotation"] =
            io::to_json(poses_t1[bodyA_id].rotation);

        std::ofstream(fmt::format("ccd-test-{:03d}.json", counter++))
            << json.dump();
    }

    void save_candidate(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1,
        const ipc::FaceVertexCandidate& fv_candidate)
    {
        nlohmann::json json;
        long bodyA_id, vertex_id, bodyB_id, face_id;
        bodies.global_to_local_vertex(
            fv_candidate.vertex_index, bodyA_id, vertex_id);
        bodies.global_to_local_face(fv_candidate.face_index, bodyB_id, face_id);
        const auto& bodyA = bodies[bodyA_id];
        const auto& bodyB = bodies[bodyB_id];
        json["type"] = "fv";

        json["face"] = nlohmann::json();
        json["face"]["vertex0"] =
            io::to_json(bodyB.vertices.row(bodyB.faces(face_id, 0)));
        json["face"]["vertex1"] =
            io::to_json(bodyB.vertices.row(bodyB.faces(face_id, 1)));
        json["face"]["vertex2"] =
            io::to_json(bodyB.vertices.row(bodyB.faces(face_id, 2)));
        json["face"]["pose_t0"] = nlohmann::json();
        json["face"]["pose_t0"]["position"] =
            io::to_json(poses_t0[bodyB_id].position);
        json["face"]["pose_t0"]["rotation"] =
            io::to_json(poses_t0[bodyB_id].rotation);
        json["face"]["pose_t1"] = nlohmann::json();
        json["face"]["pose_t1"]["position"] =
            io::to_json(poses_t1[bodyB_id].position);
        json["face"]["pose_t1"]["rotation"] =
            io::to_json(poses_t1[bodyB_id].rotation);

        json["vertex"] = nlohmann::json();
        json["vertex"]["vertex"] = io::to_json(bodyA.vertices.row(vertex_id));
        json["vertex"]["pose_t0"]["position"] =
            io::to_json(poses_t0[bodyA_id].position);
        json["vertex"]["pose_t0"]["rotation"] =
            io::to_json(poses_t0[bodyA_id].rotation);
        json["vertex"]["pose_t1"] = nlohmann::json();
        json["vertex"]["pose_t1"]["position"] =
            io::to_json(poses_t1[bodyA_id].position);
        json["vertex"]["pose_t1"]["rotation"] =
            io::to_json(poses_t1[bodyA_id].rotation);

        std::ofstream(fmt::format("ccd-test-{:03d}.json", counter++))
            << json.dump();
    }

    void save_candidate(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1,
        const ipc::EdgeEdgeCandidate& ee_candidate)
    {
        nlohmann::json json;
        long bodyA_id, edgeA_id, bodyB_id, edgeB_id;
        bodies.global_to_local_edge(
            ee_candidate.edge0_index, bodyA_id, edgeA_id);
        bodies.global_to_local_edge(
            ee_candidate.edge1_index, bodyB_id, edgeB_id);
        const auto& bodyA = bodies[bodyA_id];
        const auto& bodyB = bodies[bodyB_id];
        json["type"] = "ee";

        json["edge0"] = nlohmann::json();
        json["edge0"]["vertex0"] =
            io::to_json(bodyA.vertices.row(bodyA.edges(edgeA_id, 0)));
        json["edge0"]["vertex1"] =
            io::to_json(bodyA.vertices.row(bodyA.edges(edgeA_id, 1)));
        json["edge0"]["pose_t0"] = nlohmann::json();
        json["edge0"]["pose_t0"]["position"] =
            io::to_json(poses_t0[bodyA_id].position);
        json["edge0"]["pose_t0"]["rotation"] =
            io::to_json(poses_t0[bodyA_id].rotation);
        json["edge0"]["pose_t1"] = nlohmann::json();
        json["edge0"]["pose_t1"]["position"] =
            io::to_json(poses_t1[bodyA_id].position);
        json["edge0"]["pose_t1"]["rotation"] =
            io::to_json(poses_t1[bodyA_id].rotation);

        json["edge1"] = nlohmann::json();
        json["edge1"]["vertex0"] =
            io::to_json(bodyB.vertices.row(bodyB.edges(edgeB_id, 0)));
        json["edge1"]["vertex1"] =
            io::to_json(bodyB.vertices.row(bodyB.edges(edgeB_id, 1)));
        json["edge1"]["pose_t0"] = nlohmann::json();
        json["edge1"]["pose_t0"]["position"] =
            io::to_json(poses_t0[bodyB_id].position);
        json["edge1"]["pose_t0"]["rotation"] =
            io::to_json(poses_t0[bodyB_id].rotation);
        json["edge1"]["pose_t1"] = nlohmann::json();
        json["edge1"]["pose_t1"]["position"] =
            io::to_json(poses_t1[bodyB_id].position);
        json["edge1"]["pose_t1"]["rotation"] =
            io::to_json(poses_t1[bodyB_id].rotation);

        std::ofstream(fmt::format("ccd-test-{:03d}.json", counter++))
            << json.dump();
    }

    bool DistanceBarrierConstraint::has_active_collisions_narrow_phase(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1,
        const ipc::Candidates& candidates) const
    {
        TrajectoryType overloaded_trajectory =
            trajectory_type == TrajectoryType::PIECEWISE_LINEAR
            ? TrajectoryType::RIGID
            : trajectory_type;

        for (const auto& ev_candidate : candidates.ev_candidates) {
            double toi;
            bool are_colliding = edge_vertex_ccd(
                bodies, poses_t0, poses_t1, ev_candidate, toi,
                overloaded_trajectory);
            if (are_colliding) {
                // save_candidate(bodies, poses_t0, poses_t1, ev_candidate);
                return true;
            }
        }
        for (const auto& fv_candidate : candidates.fv_candidates) {
            double toi;
            bool are_colliding = face_vertex_ccd(
                bodies, poses_t0, poses_t1, fv_candidate, toi,
                overloaded_trajectory);
            if (are_colliding) {
                // save_candidate(bodies, poses_t0, poses_t1, fv_candidate);
                return true;
            }
        }
        for (const auto& ee_candidate : candidates.ee_candidates) {
            double toi;
            bool are_colliding = edge_edge_ccd(
                bodies, poses_t0, poses_t1, ee_candidate, toi,
                overloaded_trajectory);
            if (are_colliding) {
                // save_candidate(bodies, poses_t0, poses_t1, ee_candidate);
                return true;
            }
        }
        return false;
    }

    double DistanceBarrierConstraint::compute_earliest_toi(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1) const
    {
        PROFILE_POINT("DistanceBarrierConstraint::compute_earliest_toi");
        PROFILE_START();
        // This function will profile itself
        ipc::Candidates candidates;
        detect_collision_candidates(
            bodies, poses_t0, poses_t1, dim_to_collision_type(bodies.dim()),
            candidates, detection_method, trajectory_type,
            /*inflation_radius=*/0);

        double earliest_toi = compute_earliest_toi_narrow_phase(
            bodies, poses_t0, poses_t1, candidates);
        PROFILE_END();

        return earliest_toi;
    }

    double DistanceBarrierConstraint::compute_earliest_toi_narrow_phase(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses_t0,
        const physics::Poses<double>& poses_t1,
        const ipc::Candidates& candidates) const
    {
        NAMED_PROFILE_POINT(
            "DistanceBarrierConstraint::compute_earliest_toi_narrow_phase",
            NARROW_PHASE);
        // NOTE: These are disabled because they are not thread safe
        // NAMED_PROFILE_POINT(
        //     "DistanceBarrierConstraint::compute_earliest_toi_narrow_phase:"
        //     "edge_edge",
        //     FV_NARROW_PHASE);
        // NAMED_PROFILE_POINT(
        //     "DistanceBarrierConstraint::compute_earliest_toi_narrow_phase:"
        //     "face_vertex",
        //     EE_NARROW_PHASE);
        // NAMED_PROFILE_POINT(
        //     "DistanceBarrierConstraint::compute_earliest_toi_narrow_phase:"
        //     "edge_vertex",
        //     EV_NARROW_PHASE);

        PROFILE_START(NARROW_PHASE);

        int collision_count = 0;
        double earliest_toi = 1;
        std::mutex earliest_toi_mutex;

        const size_t num_ev = candidates.ev_candidates.size();
        const size_t num_ee = candidates.ee_candidates.size();
        const size_t num_fv = candidates.fv_candidates.size();

        // Do a single block range over all three candidate vectors
        tbb::parallel_for(
            tbb::blocked_range<int>(0, candidates.size()),
            [&](tbb::blocked_range<int> r) {
                for (int i = r.begin(); i < r.end(); i++) {
                    double toi = std::numeric_limits<double>::infinity();
                    bool are_colliding;

                    if (i < num_ev) {
                        // PROFILE_START(EV_NARROW_PHASE);
                        are_colliding = edge_vertex_ccd(
                            bodies, poses_t0, poses_t1,
                            candidates.ev_candidates[i], toi, trajectory_type,
                            earliest_toi);
                        // PROFILE_END(EV_NARROW_PHASE);
                    } else if (i - num_ev < num_ee) {
                        // PROFILE_START(EE_NARROW_PHASE);
                        are_colliding = edge_edge_ccd(
                            bodies, poses_t0, poses_t1,
                            candidates.ee_candidates[i - num_ev], toi,
                            trajectory_type, earliest_toi);
                        // PROFILE_END(EE_NARROW_PHASE);
                    } else {
                        assert(i - num_ev - num_ee < num_fv);
                        // PROFILE_START(FV_NARROW_PHASE);
                        are_colliding = face_vertex_ccd(
                            bodies, poses_t0, poses_t1,
                            candidates.fv_candidates[i - num_ev - num_ee], toi,
                            trajectory_type, earliest_toi);
                        // PROFILE_END(FV_NARROW_PHASE);
                    }

                    if (are_colliding) {
                        std::scoped_lock lock(earliest_toi_mutex);
                        collision_count++;
                        if (toi < earliest_toi) {
                            earliest_toi = toi;
                        }
                    }
                }
            });

        double percent_correct = candidates.size() == 0
            ? 100
            : (double(collision_count) / candidates.size() * 100);
        PROFILE_MESSAGE(
            NARROW_PHASE, "num_candidates,num_collisions,percentage",
            fmt::format(
                "{:d},{:d},{:g}%", candidates.size(), collision_count,
                percent_correct));

        spdlog::debug(
            "num_candidates={:d} num_collisions={:d} percentage={:g}%",
            candidates.size(), collision_count, percent_correct);

        PROFILE_END(NARROW_PHASE);

        return collision_count ? earliest_toi
                               : std::numeric_limits<double>::infinity();
    }

    void DistanceBarrierConstraint::construct_constraint_set(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses,
        ipc::Constraints& constraint_set) const
    {
        static physics::Poses<double> cached_poses;
        static ipc::Constraints cached_constraint_set;

        if (bodies.num_bodies() <= 1) {
            return;
        }

        if (poses == cached_poses) {
            constraint_set = cached_constraint_set;
            return;
        }

        PROFILE_POINT("DistanceBarrierConstraint::construct_constraint_set");
        PROFILE_START();

        const double& dhat = m_barrier_activation_distance;

        ipc::Candidates candidates;
        detect_collision_candidates_rigid(
            bodies, poses, dim_to_collision_type(bodies.dim()), candidates,
            detection_method, dhat);

        Eigen::MatrixXd V = bodies.world_vertices(poses);
        ipc::construct_constraint_set(
            candidates, /*V_rest=*/V, V, bodies.m_edges, bodies.m_faces,
            /*dhat=*/m_barrier_activation_distance, constraint_set);

        PROFILE_END();

        cached_poses = poses;
        cached_constraint_set = constraint_set;
    }

    double DistanceBarrierConstraint::compute_minimum_distance(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double>& poses) const
    {
        PROFILE_POINT("DistanceBarrierConstraint::compute_minimum_distance");
        PROFILE_START();

        ipc::Constraints constraint_set;
        construct_constraint_set(bodies, poses, constraint_set);
        Eigen::MatrixXd V = bodies.world_vertices(poses);
        double minimum_distance = sqrt(ipc::compute_minimum_distance(
            V, bodies.m_edges, bodies.m_faces, constraint_set));

        PROFILE_END();

        return minimum_distance;
    }

} // namespace opt
} // namespace ccd
