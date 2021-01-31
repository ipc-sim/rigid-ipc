#pragma once

#include <Eigen/Dense>

namespace ipc::rigid {

///
/// Generate a canonical triangle/quad subdivided from a regular grid
///
/// @param[in]  n  			 { n grid quads }
/// @param[in]  tri			 { is a tri or a quad }
/// @param[out] V            { #V x 2 output vertices positions }
/// @param[out] F            { #F x 3 output triangle indices }
///
void regular_2d_grid(
    const int n, const bool tri, Eigen::MatrixXd& V, Eigen::MatrixXi& F);

void regular_2d_grid(
    const int num_cols,
    const int num_rows,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F);

} // namespace ipc::rigid
