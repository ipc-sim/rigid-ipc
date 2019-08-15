#include "collision_constraint.hpp"

#include <profiler.hpp>

namespace ccd {
// !important: this needs to be define in the enum namespace
NLOHMANN_JSON_SERIALIZE_ENUM(DetectionMethod,
    { { HASH_GRID, "hash_grid" }, { BRUTE_FORCE, "brute_force" } })

namespace opt {

    CollisionConstraint::CollisionConstraint(const std::string& name)
        : detection_method(HASH_GRID)
        , name_(name)
    {
    }

    CollisionConstraint::~CollisionConstraint() {}

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

    void CollisionConstraint::initialize(const Eigen::MatrixX2d& V,
        const Eigen::MatrixX2i& E,
        const Eigen::VectorXi& Gid,
        const Eigen::MatrixXd& Uk)
    {
        vertices = V;
        edges = E;
        group_ids = Gid;

        m_ev_impacts.clear();
        update_collision_set(Uk);
    }

    void CollisionConstraint::update_collision_set(const Eigen::MatrixXd& Uk)
    {
        ccd::detect_edge_vertex_collisions(vertices, Uk, edges, group_ids,
            m_ev_impacts, detection_method);
    }

} // namespace opt
} // namespace ccd
