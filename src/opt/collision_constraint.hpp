#pragma once

#include <nlohmann/json.hpp>

#include <ccd/rigid_body_collision_detection.hpp>
#include <physics/rigid_body_assembler.hpp>

namespace ccd {
namespace opt {

    class CollisionConstraint {
    public:
        CollisionConstraint(const std::string& name);
        virtual ~CollisionConstraint() = default;

        virtual void settings(const nlohmann::json& json);
        virtual nlohmann::json settings() const;

        inline const std::string& name() const { return name_; }

        virtual void initialize() {};

        void construct_collision_set(
            const physics::RigidBodyAssembler& bodies,
            const physics::Poses<double> poses_t0,
            const physics::Poses<double> poses_t1,
            ConcurrentImpacts& impacts) const;

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

        std::string name_;
    };

} // namespace opt
} // namespace ccd
