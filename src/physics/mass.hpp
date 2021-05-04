#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>

#include <utils/eigen_ext.hpp>

namespace ipc::rigid {

/// @brief Compute the total mass, center of mass, and moment of intertia
void compute_mass_properties(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& facets,
    double& total_mass,
    VectorMax3d& center_of_mass,
    MatrixMax3d& moment_of_inertia);

/// @brief Compute the 2D total mass, center of mass, and moment of intertia
void compute_mass_properties_2D(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& edges,
    double& mass,
    VectorMax3d& center,
    MatrixMax3d& intertia);

/// @brief Compute the 3D total mass, center of mass, and moment of intertia
void compute_mass_properties_3D(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& faces,
    double& mass,
    VectorMax3d& center,
    MatrixMax3d& intertia);

/// @brief Construct the sparse mass matrix for the given mesh (V, E).
void construct_mass_matrix(
    const Eigen::MatrixXd& vertices,
    const Eigen::MatrixXi& facets,
    Eigen::SparseMatrix<double>& mass_matrix);

/// @brief Computes the total mass for the given mesh
double compute_total_mass(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& facets);
/// @brief Computes the total mass from the mass matrix
double compute_total_mass(const Eigen::SparseMatrix<double>& mass_matrix);

VectorMax3d compute_center_of_mass(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& facets);
VectorMax3d compute_center_of_mass(
    const Eigen::MatrixXd& vertices,
    const Eigen::SparseMatrix<double>& mass_matrix);

/**
 * @brief Computes the moment of intertia
 *
 * Assumes vertices are given in body space (i.e centered of mass at 0,0).
 */
MatrixMax3d compute_moment_of_inertia(
    const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& facets);
MatrixMax3d compute_moment_of_inertia(
    const Eigen::MatrixXd& vertices,
    const Eigen::SparseMatrix<double>& mass_matrix);

} // namespace ipc::rigid
