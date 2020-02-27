#include "collision_constraint.hpp"

namespace ccd {
// !important: this needs to be define in the enum namespace
NLOHMANN_JSON_SERIALIZE_ENUM(
    DetectionMethod,
    { { HASH_GRID, "hash_grid" }, { BRUTE_FORCE, "brute_force" } })

namespace opt {

    CollisionConstraint::CollisionConstraint(const std::string& name)
        : detection_method(HASH_GRID)
        , name_(name)
    {
    }

    void CollisionConstraint::settings(const nlohmann::json& json)
    {
        detection_method = json["detection_method"].get<DetectionMethod>();
    }

    nlohmann::json CollisionConstraint::settings() const
    {
        nlohmann::json json;
        json["detection_method"] = detection_method;
        return json;
    }

    void CollisionConstraint::construct_collision_set(
        const physics::RigidBodyAssembler& bodies,
        const physics::Poses<double> poses,
        const physics::Poses<double> displacements,
        ConcurrentImpacts& impacts) const
    {
        ccd::detect_collisions(
            bodies, poses, displacements, dim_to_collision_type(bodies.dim()),
            impacts, detection_method);
    }

} // namespace opt
} // namespace ccd
