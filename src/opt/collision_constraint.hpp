#pragma once

#include <nlohmann/json.hpp>

#include <ccd/ccd.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ipc::rigid {

class CollisionConstraint {
public:
    CollisionConstraint(const std::string& name);
    virtual ~CollisionConstraint() = default;

    virtual void settings(const nlohmann::json& json);
    virtual nlohmann::json settings() const;

    inline const std::string& name() const { return m_name; }

    virtual void initialize() {};

    void construct_collision_set(
        const RigidBodyAssembler& bodies,
        const PosesD poses_t0,
        const PosesD poses_t1,
        Impacts& impacts) const;

    // Settings
    // ----------
    DetectionMethod detection_method;
    TrajectoryType trajectory_type;

protected:
    inline static int dim_to_collision_type(int dim)
    {
        return dim == 2
            ? CollisionType::EDGE_VERTEX
            : (CollisionType::EDGE_EDGE | CollisionType::FACE_VERTEX);
    }

    std::string m_name;
};

} // namespace ipc::rigid
