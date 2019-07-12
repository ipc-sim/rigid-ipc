#pragma once

#include <igl/opengl/ViewerData.h>
#include <igl/opengl/glfw/Viewer.h>

namespace igl {
namespace opengl {

    class ViewerDataExt {
    public:
        ViewerDataExt(igl::opengl::glfw::Viewer* _viewer);
        ViewerData& data() { return m_viewer->data_list[data_id]; }

        void set_vector_field(const Eigen::MatrixXd& V,
            const Eigen::MatrixXd& F,
            const Eigen::RowVector3d& color);

        void set_graph(const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& E,
            const Eigen::RowVector3d& color);

        void set_vertex_data(
            const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> data);
        void update_vertex_data();
        void update_graph(const Eigen::MatrixXd& V);
        void update_vector_field(
            const Eigen::MatrixXd& V, const Eigen::MatrixXd& F);

        void recolor();
        inline bool is_graph() { return mE.size() > 0; }


        igl::opengl::glfw::Viewer* m_viewer;
        size_t data_id;

        Eigen::MatrixXi mE;
        Eigen::MatrixXd mV;
        Eigen::MatrixXd mF;
        Eigen::RowVector3d m_color;
        std::vector<std::string> vertex_data_labels;
        bool show_vertex_data;

    protected:
        Eigen::MatrixXd set_vertices(const Eigen::MatrixXd& V);

        void set_points(const Eigen::MatrixXd& V, const Eigen::MatrixXd& color);

        void set_edges(const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& E,
            const Eigen::MatrixXd& color);

        void set_vector_field(const Eigen::MatrixXd& V,
            const Eigen::MatrixXd& F,
            const Eigen::MatrixXd& color);
    };

} // namespace opengl
} // namespace igl
