#include "center_of_mass.hpp"

namespace ccd {
namespace physics {

    Eigen::VectorXd center_of_mass(
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges)
    {
        assert(edges.cols() == 2);
        assert(vertices.cols() == 2);

        Eigen::VectorXd vertex_masses = Eigen::VectorXd::Zero(vertices.rows());
        double total_edge_length = 0;
        for (long i = 0; i < edges.rows(); i++) {
            double edge_length
                = (vertices.row(edges(i, 1)) - vertices.row(edges(i, 0)))
                      .norm();
            vertex_masses(edges(i, 0)) += edge_length / 2;
            vertex_masses(edges(i, 1)) += edge_length / 2;
            total_edge_length += edge_length;
        }
        Eigen::VectorXd cm;
        cm = (vertex_masses.asDiagonal() * vertices).colwise().sum()
            / total_edge_length;

        return cm;
    }
} // namespace physics
} // namespace ccd
