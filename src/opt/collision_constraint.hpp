#pragma once

#include <Eigen/SparseCore>
#include <memory> // shared_ptr

#include <nlohmann/json.hpp>

#include <ccd/collision_detection.hpp>
#include <utils/not_implemented_error.hpp>

#include <autodiff/autodiff_types.hpp>

namespace ccd {
namespace opt {


    class CollisionConstraint {
    public:
        CollisionConstraint(const std::string& name);
        virtual ~CollisionConstraint();

        virtual void settings(const nlohmann::json& json);
        virtual nlohmann::json settings() const;

        inline const std::string& name() const { return name_; }

        virtual void initialize(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges,
            const Eigen::VectorXi& group_ids,
            const Eigen::MatrixXd& Uk);

        virtual void update_collision_set(const Eigen::MatrixXd& Uk);
        const EdgeVertexImpacts& ev_impacts() {return m_ev_impacts; }

        // Settings
        // ----------
        DetectionMethod detection_method;

        // Structures used for detection
        // ------------
        Eigen::MatrixX2d vertices;
        Eigen::MatrixX2i edges;
        Eigen::VectorXi group_ids;

    protected:
        /// @brief All edge-vertex contact
        EdgeVertexImpacts m_ev_impacts;
        std::string name_;
    };

    void slice_vector(const Eigen::MatrixX2d& data,
        const Eigen::Vector2i e_ij,
        Eigen::Vector2i e_kl,
        std::array<Eigen::Vector2d, 4>& d);

    template <typename T> bool differentiable();

} // namespace opt
} // namespace ccd
