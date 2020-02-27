#pragma once

#include <nlohmann/json.hpp>

#include <ccd/collision_detection.hpp>

namespace ccd {
namespace opt {


    class CollisionConstraint {
    public:
        CollisionConstraint(const std::string& name);
        virtual ~CollisionConstraint() = default;

        virtual void settings(const nlohmann::json& json);
        virtual nlohmann::json settings() const;

        inline const std::string& name() const { return name_; }

        virtual EdgeVertexImpacts initialize(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges,
            const Eigen::VectorXi& group_ids,
            const Eigen::MatrixXd& Uk);

        EdgeVertexImpacts get_collision_set(const Eigen::MatrixXd& Uk);

        // Settings
        // ----------
        DetectionMethod detection_method;

        // Structures used for detection
        // ------------
        Eigen::MatrixX2d vertices;
        Eigen::MatrixX2i edges;
        Eigen::VectorXi group_ids;

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
