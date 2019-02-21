#include "edges_to_rectangles.hpp"

namespace ccd {
namespace widgets {

    void edges_to_rectangles(
        const Eigen::MatrixX2d& vertex,
        const Eigen::MatrixX2i& edges,
        const Eigen::VectorXd& width,
        const double min_width,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        const bool nodes_only)
    {
        long num_edges = edges.rows();
        Eigen::MatrixX2d nodes_s(num_edges, 2), nodes_f(num_edges, 2);

        for (int i = 0; i < num_edges; ++i) {
            nodes_s.row(i) << vertex.row(edges(i, 0));
            nodes_f.row(i) << vertex.row(edges(i, 1));
        }

        return lines_to_rectangles(nodes_s, nodes_f, width, min_width, V, F, nodes_only);
    }

    void lines_to_rectangles(
        const Eigen::MatrixX2d& nodes_s,
        const Eigen::MatrixX2d& nodes_f,
        const Eigen::VectorXd& width,
        const double min_width,
        Eigen::MatrixXd& V,
        Eigen::MatrixXi& F,
        const bool nodes_only)
    {
        assert(nodes_s.rows() == nodes_f.rows());
        assert(nodes_s.rows() == width.rows());

        int num_rect = int(nodes_s.rows());
        int num_nodes = num_rect * 4;
        int num_tris = num_rect * 2;

        V.resize(num_nodes, 3);
        V.setZero();

        Eigen::MatrixX2d e = nodes_f - nodes_s;
        Eigen::ArrayXd e_len = e.rowwise().norm();
        Eigen::ArrayX2d e_normal(num_rect, 2);
        e_normal.col(0) = -e.col(1);
        e_normal.col(1) = e.col(0);
        e_normal = e_normal.colwise() / e_len;

        Eigen::ArrayX2d n0 = nodes_s.array();
        Eigen::ArrayX2d n1 = nodes_f.array();
        Eigen::ArrayXd w = width.array() / 2.0;
        w = (w < min_width).select(min_width, w);

        V.block(0 * num_rect, 0, num_rect, 2) = n0 + e_normal.colwise() * w;
        V.block(1 * num_rect, 0, num_rect, 2) = n0 - e_normal.colwise() * w;
        V.block(2 * num_rect, 0, num_rect, 2) = n1 - e_normal.colwise() * w;
        V.block(3 * num_rect, 0, num_rect, 2) = n1 + e_normal.colwise() * w;

        if (nodes_only){
            return;
        }

        // face 1: 0, 1, 2
        // face 2: 2, 3, 0
        Eigen::MatrixX3i F0(num_rect, 3), F1(num_rect, 3);
        Eigen::RowVector3i f0, f1;
        f0 << 0 * num_rect, 1 * num_rect, 2 * num_rect;
        f1 << 2 * num_rect, 3 * num_rect, 0 * num_rect;

        F0.rowwise() = f0;
        F1.rowwise() = f1;

        Eigen::ArrayXi x = Eigen::ArrayXi::LinSpaced(num_rect, 0, num_rect - 1);

        F.resize(num_tris, 3);
        F.block(0 * num_rect, 0, num_rect, 3) = F0.array().colwise() + x;
        F.block(1 * num_rect, 0, num_rect, 3) = F1.array().colwise() + x;
    }

}
}
