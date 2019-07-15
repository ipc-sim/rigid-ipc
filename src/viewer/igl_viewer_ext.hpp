#pragma once

#include <igl/opengl/ViewerData.h>
#include <igl/opengl/glfw/Viewer.h>

#include <igl/colormap.h>

namespace igl {
namespace opengl {

    class ViewerDataExt {
    public:
        ViewerDataExt(igl::opengl::glfw::Viewer* _viewer);
        virtual ~ViewerDataExt();
        ViewerData& data() { return m_viewer->data_list[data_id]; }

        virtual void recolor() {}
        virtual inline bool is_graph() { return false; }
        virtual inline bool is_scalar_field() { return false; }

        /// only used in graph data for now
        virtual void update_vertex_data() {}

        igl::opengl::glfw::Viewer* m_viewer;
        size_t data_id;

        Eigen::RowVector3d m_color;

        bool show_vertex_data;
        bool m_log_scale;

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

    class GraphData : public ViewerDataExt {
    public:
        GraphData(igl::opengl::glfw::Viewer* _viewer);

        void set_graph(const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& E,
            const Eigen::RowVector3d& color);

        void set_vertex_data(
            const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> data);
        void update_vertex_data() override;
        void update_graph(const Eigen::MatrixXd& V);

        void recolor() override;
        inline bool is_graph() override { return true; }
        Eigen::MatrixXi mE;
        Eigen::MatrixXd mV;

        std::vector<std::string> vertex_data_labels;
    };

    class VectorFieldData : public ViewerDataExt {
    public:
        VectorFieldData(igl::opengl::glfw::Viewer* _viewer);

        void recolor() override;
        void set_vector_field(const Eigen::MatrixXd& V,
            const Eigen::MatrixXd& F,
            const Eigen::RowVector3d& color);

        void update_vector_field(
            const Eigen::MatrixXd& V, const Eigen::MatrixXd& F);

        Eigen::MatrixXd mF;
        Eigen::MatrixXd mV;
    };

    class ScalarFieldData : public ViewerDataExt {
    public:
        ScalarFieldData(igl::opengl::glfw::Viewer* _viewer);

        void recolor() override;
        void set_mesh(const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& F,
            const Eigen::RowVector3d& base_color,
            const Eigen::RowVector3d& inf_color);

        void set_vertex_data(const Eigen::MatrixXd& vertex_data);

        inline bool is_scalar_field() override { return true; }

        Eigen::MatrixXi mF;
        Eigen::MatrixXd mV;
        Eigen::MatrixXd m_vertex_data;

        igl::ColorMapType m_colormap_type;
        Eigen::Vector3d m_base_color;

    };

} // namespace opengl
} // namespace igl
