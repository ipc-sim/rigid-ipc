#include "rigid_body.hpp"

#include <Eigen/Eigenvalues>
#include <Eigen/Geometry>

#include <autodiff/autodiff_types.hpp>
#include <finitediff.hpp>
#include <logger.hpp>
#include <physics/mass.hpp>
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
        const int group_id)
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

        double m;
        Eigen::VectorX3d center_of_mass;
        Eigen::MatrixXX3d I;
        compute_mass_properties(
            vertices_, dim == 2 || faces.size() == 0 ? edges : faces, m,
            center_of_mass, I);
        Eigen::MatrixXd centered_vertices =
            vertices_.rowwise() - center_of_mass.transpose();

        // set position so current vertices match input
        Pose<double> adjusted_pose(center_of_mass, pose.rotation);

        assert(is_dof_fixed.size() == pose.ndof());
        return RigidBody(
            centered_vertices, edges, faces, adjusted_pose, velocity, force,
            density, is_dof_fixed, oriented, group_id);
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
        const int group_id)
        : group_id(group_id)
        , vertices(vertices)
        , edges(edges)
        , faces(faces)
        , is_dof_fixed(is_dof_fixed)
        , is_oriented(oriented)
        , pose(pose)
        , pose_prev(pose)
        , velocity(velocity)
        , velocity_prev(velocity)
        , force(force)
    {
        assert(edges.size() == 0 || edges.cols() == 2);
        assert(faces.size() == 0 || faces.cols() == 3);

        Eigen::VectorX3d center_of_mass;
        Eigen::MatrixXX3d I;
        compute_mass_properties(
            vertices, dim() == 2 || faces.size() == 0 ? edges : faces, mass,
            center_of_mass, I);
        assert(center_of_mass.squaredNorm() < 1e-8);

        // TODO: Not sure why this is times based on Chrono
        // (https://bit.ly/2TVjJVm). Might be because mass above is actually
        // volume.
        mass *= density;
        if (dim() == 3) {
            // Got this from Chrono: https://bit.ly/2RpbTl1
            Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> es;
            es.compute(I);
            moment_of_inertia = density * es.eigenvalues();
            R0 = es.eigenvectors();
            // Ensure that we have an orientation preserving transform
            if (R0.determinant() < 0.0) {
                R0.col(0) *= -1.0;
            }
            assert(R0.isUnitary(1e-9));
            assert(fabs(R0.determinant() - 1.0) <= 1.0e-9);
            Eigen::VectorX3<bool> is_rot_dof_fixed =
                is_dof_fixed.tail(Pose<double>::dim_to_rot_ndof(dim()));
            if (is_rot_dof_fixed.count() == 2) {
                // Convert moment of inertia to world coordinates
                // https://physics.stackexchange.com/a/268812
                moment_of_inertia = -I.diagonal().array() + I.diagonal().sum();
                R0.setIdentity();
            } else if (is_rot_dof_fixed.count() == 1) {
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
            // τ = R₀ᵀτ₀ (τ₀ expressed in body coordinates)
            this->force.rotation = R0.transpose() * this->force.rotation;
        } else {
            moment_of_inertia = density * I.diagonal();
            R0 = Eigen::Matrix<double, 1, 1>::Identity();
        }

        // Zero out the velocity and forces of fixed dof
        this->velocity.zero_dof(is_dof_fixed, R0);
        this->force.zero_dof(is_dof_fixed, R0);

        // Compute and construct some useful constants
        mass_matrix.resize(ndof());
        mass_matrix.diagonal().head(pos_ndof()).setConstant(mass);
        mass_matrix.diagonal().tail(rot_ndof()) = moment_of_inertia;

        r_max = vertices.rowwise().squaredNorm().maxCoeff();

        average_edge_length = 0;
        for (long i = 0; i < edges.rows(); i++) {
            average_edge_length +=
                (vertices.row(edges(i, 0)) - vertices.row(edges(i, 1))).norm();
        }
        average_edge_length /= edges.rows();
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

    Eigen::MatrixXd
    RigidBody::world_vertices_gradient(const Pose<double>& _pose) const
    {
        // Dynamic number of dof, but there is a limit of 6 dof in 3D.
        typedef AutodiffType<Eigen::Dynamic, 6> Diff;
        Diff::activate(_pose.ndof());

        Pose<Diff::DDouble1> dpose(
            Diff::d1vars(0, _pose.position),
            Diff::d1vars(_pose.pos_ndof(), _pose.rotation));
        dpose.rotation = dpose.rotation / Diff::DDouble1(r_max);

        Diff::D1MatrixXd dx = world_vertices<Diff::DDouble1>(dpose);

        flatten<Diff::DDouble1>(dx);
        Eigen::MatrixXd gradient = Diff::get_gradient(dx);
        return gradient;
    }

} // namespace physics
} // namespace ccd
