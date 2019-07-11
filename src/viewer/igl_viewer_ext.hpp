#pragma once

#include <igl/opengl/ViewerData.h>
namespace igl {
namespace opengl {

    class ViewerDataExt {
    public:
        ViewerDataExt(ViewerData& data);

        void set_graph(const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& E,
            const Eigen::RowVector3d& color);

        void set_vertex_data(const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> data);
        void update_graph(const Eigen::MatrixXd& V);

        ViewerData& m_data;
        Eigen::MatrixXi mE;
        Eigen::MatrixXd mV;
        Eigen::RowVector3d m_color;
        std::vector<std::string> vertex_data_labels;

    protected:
        Eigen::MatrixXd set_vertices(const Eigen::MatrixXd& V);

        void set_points(const Eigen::MatrixXd& V, const Eigen::MatrixXd& color);

        void set_edges(const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& E,
            const Eigen::MatrixXd& color);
    };

} // namespace opengl
} // namespace igl
