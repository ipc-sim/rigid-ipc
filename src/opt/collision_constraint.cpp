#include "collision_constraint.hpp"

namespace ipc::rigid {

CollisionConstraint::CollisionConstraint(const std::string& name)
    : detection_method(DetectionMethod::HASH_GRID)
    , trajectory_type(TrajectoryType::RIGID)
    , m_name(name)
{
}

void CollisionConstraint::settings(const nlohmann::json& json)
{
    detection_method = json["detection_method"].get<DetectionMethod>();
    trajectory_type = json["trajectory_type"].get<TrajectoryType>();
}

nlohmann::json CollisionConstraint::settings() const
{
    nlohmann::json json;
    json["detection_method"] = detection_method;
    json["trajectory_type"] = trajectory_type;
    return json;
}

void CollisionConstraint::construct_collision_set(
    const RigidBodyAssembler& bodies,
    const PosesD poses_t0,
    const PosesD poses_t1,
    Impacts& impacts) const
{
    detect_collisions(
        bodies, poses_t0, poses_t1, dim_to_collision_type(bodies.dim()),
        impacts, detection_method, trajectory_type);
}

} // namespace ipc::rigid
