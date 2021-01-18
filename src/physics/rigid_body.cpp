#include "rigid_body.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <autodiff/autodiff_types.hpp>
#include <finitediff.hpp>
#include <logger.hpp>
#include <physics/mass.hpp>
#include <profiler.hpp>
#include <utils/eigen_ext.hpp>
#include <utils/flatten.hpp>
#include <utils/not_implemented_error.hpp>

namespace ccd {

namespace physics {

    RigidBody RigidBody::from_points(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const Pose<double>& pose,
        const Pose<double>& velocity,
        const Pose<double>& force,
        const double density,
        const Eigen::VectorX6b& is_dof_fixed,
        const bool oriented,
        const int group_id,
        const RigidBodyType type,
        const double kinematic_max_time)
    {
        int dim = vertices.cols();
        assert(dim == pose.dim());
        assert(dim == velocity.dim());
        assert(dim == force.dim());
        assert(edges.size() == 0 || edges.cols() == 2);
        assert(faces.size() == 0 || faces.cols() == 3);

        // move vertices so their center of mass is at (0, 0)
        Eigen::MatrixXd vertices_ = vertices;
        vertices_.rowwise() += pose.position.transpose();

        // compute the center of mass several times to get more accurate
        Eigen::MatrixXd centered_vertices = vertices_;
        Eigen::VectorX3d center_of_mass = Eigen::VectorX3d::Zero(dim);
        for (int i = 0; i < 10; i++) {
            double tmp_m;
            Eigen::VectorX3d tmp_center_of_mass;
            Eigen::MatrixXX3d tmp_I;
            compute_mass_properties(
                centered_vertices,
                dim == 2 || faces.size() == 0 ? edges : faces, //
                tmp_m, tmp_center_of_mass, tmp_I);
            centered_vertices.rowwise() -= tmp_center_of_mass.transpose();
            center_of_mass += tmp_center_of_mass;
            if (tmp_center_of_mass.squaredNorm() < 1e-8) {
                break;
            }
        }

        // set position so current vertices match input
        Pose<double> adjusted_pose(center_of_mass, pose.rotation);

        assert(is_dof_fixed.size() == pose.ndof());
        return RigidBody(
            centered_vertices, edges, faces, adjusted_pose, velocity, force,
            density, is_dof_fixed, oriented, group_id, type,
            kinematic_max_time);
    }

    RigidBody::RigidBody(
        const Eigen::MatrixXd& vertices,
        const Eigen::MatrixXi& edges,
        const Eigen::MatrixXi& faces,
        const Pose<double>& pose,
        const Pose<double>& velocity,
        const Pose<double>& force,
        const double density,
        const Eigen::VectorX6b& is_dof_fixed,
        const bool oriented,
        const int group_id,
        const RigidBodyType type,
        const double kinematic_max_time)
        : group_id(group_id)
        , type(type)
        , kinematic_max_time(kinematic_max_time)
        , vertices(vertices)
        , edges(edges)
        , faces(faces)
        , is_dof_fixed(is_dof_fixed)
        , is_oriented(oriented)
        , mesh_selector(vertices.rows(), edges, faces)
        , pose(pose)
        , velocity(velocity)
        , force(force)
    {
        assert(edges.size() == 0 || edges.cols() == 2);
        assert(faces.size() == 0 || faces.cols() == 3);

        if (type == RigidBodyType::STATIC) {
            this->is_dof_fixed.setOnes(this->is_dof_fixed.size());
        } else if (this->is_dof_fixed.array().all()) {
            this->type = RigidBodyType::STATIC;
        }

        Eigen::VectorX3d center_of_mass;
        Eigen::MatrixXX3d I;
        compute_mass_properties(
            vertices, dim() == 2 || faces.size() == 0 ? edges : faces, mass,
            center_of_mass, I);
        // assert(center_of_mass.squaredNorm() < 1e-8);

        // Mass above is actually volume in m³ and density is Kg/m³
        mass *= density;
        if (dim() == 3) {
            // Got this from Chrono: https://bit.ly/2RpbTl1
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
            double threshold = I.lpNorm<Eigen::Infinity>() * 1e-16;
            I = (threshold < I.array().abs()).select(I, 0.0);
            es.compute(I);
            if (es.info() != Eigen::Success) {
                spdlog::error(
                    "Eigen decompostion of the inertia tensor failed!");
            }
            moment_of_inertia = density * es.eigenvalues();
            if ((moment_of_inertia.array() < 0).any()) {
                spdlog::warn(
                    "Negative moment of inertia ({}), projecting to 0.",
                    logger::fmt_eigen(moment_of_inertia));
                // Avoid negative epsilon inertias
                moment_of_inertia = (moment_of_inertia.array() < 0)
                                        .select(0, moment_of_inertia);
            }
            R0 = es.eigenvectors();
            // Ensure that we have an orientation preserving transform
            if (R0.determinant() < 0.0) {
                R0.col(0) *= -1.0;
            }
            assert(R0.isUnitary(1e-9));
            assert(fabs(R0.determinant() - 1.0) <= 1.0e-9);
            int num_rot_dof_fixed =
                is_dof_fixed.tail(Pose<double>::dim_to_rot_ndof(dim())).count();
            if (num_rot_dof_fixed == 2) {
                // Convert moment of inertia to world coordinates
                // https://physics.stackexchange.com/a/268812
                moment_of_inertia = -I.diagonal().array() + I.diagonal().sum();
                R0.setIdentity();
            } else if (num_rot_dof_fixed == 1) {
                spdlog::warn("Rigid body dynamics with two rotational DoF has "
                             "not been tested thoroughly.");
            }
            // R = RᵢR₀
            Eigen::AngleAxisd r = Eigen::AngleAxisd(
                Eigen::Matrix3d(this->pose.construct_rotation_matrix() * R0));
            this->pose.rotation = r.angle() * r.axis();
            // v = Rv₀ + p = RᵢR₀v₀ + p = RᵢR₀R₀ᵀv₀ + p
            this->vertices = this->vertices * R0; // R₀ᵀ * V₀ᵀ = V₀ * R₀
            // ω = R₀ᵀω₀ (ω₀ expressed in body coordinates)
            this->velocity.rotation = R0.transpose() * this->velocity.rotation;
            Eigen::Matrix3d Q_t0 = this->pose.construct_rotation_matrix();
            this->Qdot = Q_t0 * Eigen::Hat(this->velocity.rotation);
            // τ = R₀ᵀτ₀ (τ₀ expressed in body coordinates)
            // NOTE: this transformation will be done later
            // this->force.rotation = R0.transpose() * this->force.rotation;
        } else {
            moment_of_inertia = density * I.diagonal();
            R0 = Eigen::Matrix<double, 1, 1>::Identity();
        }

        // Zero out the velocity and forces of fixed dof
        this->velocity.zero_dof(is_dof_fixed, R0);
        this->force.zero_dof(is_dof_fixed, R0);

        // Update the previous pose and velocity to reflect the changes made
        // here
        this->pose_prev = this->pose;
        this->velocity_prev = this->velocity;

        // Compute and construct some useful constants
        mass_matrix.resize(ndof());
        mass_matrix.diagonal().head(pos_ndof()).setConstant(mass);
        mass_matrix.diagonal().tail(rot_ndof()) = moment_of_inertia;

        r_max = vertices.rowwise().norm().maxCoeff();

        average_edge_length = 0;
        for (long i = 0; i < edges.rows(); i++) {
            average_edge_length +=
                (vertices.row(edges(i, 0)) - vertices.row(edges(i, 1))).norm();
        }
        if (edges.rows() > 0) {
            average_edge_length /= edges.rows();
        }
        assert(std::isfinite(average_edge_length));

        init_bvh();
    }

    void RigidBody::init_bvh()
    {
        PROFILE_POINT("RigidBody::init_bvh");
        PROFILE_START();

        // heterogenous bounding boxes
        std::vector<std::array<Eigen::Vector3d, 2>> aabbs(
            num_codim_vertices() + num_codim_edges() + num_faces());

        for (size_t i = 0; i < num_codim_vertices(); i++) {
            size_t vi = mesh_selector.codim_vertices_to_vertices(i);
            if (dim() == 2) {
                aabbs[i][0][2] = 0;
                aabbs[i][1][2] = 0;
            }
            aabbs[i][0].head(dim()) = vertices.row(i);
            aabbs[i][1].head(dim()) = vertices.row(i);
        }

        size_t start_i = num_codim_vertices();
        for (size_t i = 0; i < num_codim_edges(); i++) {
            size_t ei = mesh_selector.codim_edges_to_edges(i);
            const auto& e0 = vertices.row(edges(ei, 0));
            const auto& e1 = vertices.row(edges(ei, 1));

            if (dim() == 2) {
                aabbs[start_i + i][0][2] = 0;
                aabbs[start_i + i][1][2] = 0;
            }
            aabbs[start_i + i][0].head(dim()) = e0.cwiseMin(e1);
            aabbs[start_i + i][1].head(dim()) = e0.cwiseMax(e1);
        }

        start_i += num_codim_edges();
        for (size_t i = 0; i < num_faces(); i++) {
            assert(dim() == 3);
            const auto& f0 = vertices.row(faces(i, 0));
            const auto& f1 = vertices.row(faces(i, 1));
            const auto& f2 = vertices.row(faces(i, 2));
            aabbs[start_i + i][0] = f0.cwiseMin(f1).cwiseMin(f2);
            aabbs[start_i + i][1] = f0.cwiseMax(f1).cwiseMax(f2);
        }

        bvh.init(aabbs);

        PROFILE_END();
    }

    Eigen::MatrixXd RigidBody::world_velocities() const
    {
        // compute ẋ = Q̇ * x_B + q̇
        // where Q̇ = Qω̂
        Eigen::MatrixXX3d dQ_dt =
            pose.construct_rotation_matrix() * Eigen::Hat(velocity.rotation);
        return (vertices * dQ_dt.transpose()).rowwise()
            + velocity.position.transpose();
    }

    void RigidBody::compute_bounding_box(
        const physics::Pose<double>& pose_t0,
        const physics::Pose<double>& pose_t1,
        Eigen::VectorX3d& box_min,
        Eigen::VectorX3d& box_max) const
    {
        PROFILE_POINT("RigidBody::compute_bounding_box");
        PROFILE_START();

        // If the body is not rotating then just use the linearized
        // trajectory
        if (type == RigidBodyType::STATIC
            || (pose_t0.rotation.array() == pose_t1.rotation.array()).all()) {
            Eigen::MatrixXd V0 = world_vertices(pose_t0);
            box_min = V0.colwise().minCoeff();
            box_max = V0.colwise().maxCoeff();
            Eigen::MatrixXd V1 = world_vertices(pose_t1);
            box_min = box_min.cwiseMin(V1.colwise().minCoeff().transpose());
            box_max = box_max.cwiseMax(V1.colwise().maxCoeff().transpose());
        } else {
            // Use the maximum radius of the body to bound all rotations
            box_min =
                pose_t0.position.cwiseMin(pose_t1.position).array() - r_max;
            box_max =
                pose_t0.position.cwiseMax(pose_t1.position).array() + r_max;
        }

        PROFILE_END();
    }

} // namespace physics
} // namespace ccd
