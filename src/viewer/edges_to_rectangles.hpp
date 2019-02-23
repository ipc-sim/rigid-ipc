#ifndef EDGES_TO_RECTANGLES_H
#define EDGES_TO_RECTANGLES_H

#include <Eigen/Core>

namespace ccd {
namespace widgets {

    void edges_to_rectangles(
        const Eigen::MatrixX2d& vertex,
        const Eigen::MatrixX2i& edges,
        const Eigen::VectorXd& width,
        const double min_width,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        bool nodes_only = false);

    void lines_to_rectangles(
        const Eigen::MatrixX2d& nodes_s,
        const Eigen::MatrixX2d& nodes_f,
        const Eigen::VectorXd& width,
        const double min_width,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        bool nodes_only = false);

}
}
#endif // EDGES_TO_RECTANGLES_H
