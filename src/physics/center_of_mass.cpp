#include "center_of_mass.hpp"

#include <physics/mass_matrix.hpp>

namespace ccd {
namespace physics {

    Eigen::VectorXd center_of_mass(
        const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& edges)
    {
        assert(edges.cols() == 2);
        assert(vertices.cols() == 2);

        Eigen::VectorXd vertex_masses;
        mass_vector(vertices,edges, vertex_masses);

        Eigen::VectorXd cm;
        cm = (vertex_masses.asDiagonal() * vertices).colwise().sum()
            / vertex_masses.sum();

        return cm;
    }
} // namespace physics
} // namespace ccd
