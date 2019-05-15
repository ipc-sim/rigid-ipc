#pragma once

#include <Eigen/SparseCore>

#include <ccd/collision_detection.hpp>

namespace ccd {
namespace opt {

    class CollisionConstraint {
    public:
        CollisionConstraint();
        virtual ~CollisionConstraint();

        virtual void initialize(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk)
            = 0;

        void detecteCollisions(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk);

        virtual void compute_constraints(const Eigen::MatrixX2d& vertices,
            const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& Uk,
            Eigen::VectorXd& g_uk)
            = 0;

        virtual void compute_constraints_jacobian(
            const Eigen::MatrixX2d& vertices, const Eigen::MatrixX2i& edges,
            const Eigen::MatrixXd& Uk, Eigen::MatrixXd& g_uk_jacobian)
            = 0;

        virtual void compute_constraints_hessian(
            const Eigen::MatrixX2d& vertices, const Eigen::MatrixX2i& edges,
            const Eigen::MatrixXd& Uk,
            std::vector<Eigen::SparseMatrix<double>>& g_uk_hessian)
            = 0;

        // Settings
        // ----------
        DetectionMethod detection_method;
        bool recompute_collision_set;

        // Structures used for detection
        // ------------
        /// @brief All edge-vertex contact
        EdgeVertexImpacts ev_impacts;

        /// @brief All edge-edge contact
        EdgeEdgeImpacts ee_impacts;

        /// @brief #E,1 indices of the edges' first impact
        Eigen::VectorXi edge_impact_map;

        /// @brief The current number of pruned impacts
        int num_pruned_impacts;
    };

} // namespace opt
} // namespace ccd
