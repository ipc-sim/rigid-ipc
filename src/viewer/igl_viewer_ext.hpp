#pragma once

#include <igl/colormap.h>
#include <igl/opengl/ViewerData.h>
#include <igl/opengl/glfw/Viewer.h>

#include <physics/pose.hpp>
#include <utils/eigen_ext.hpp>

using namespace ipc;

namespace igl {
namespace opengl {

    class ViewerDataExt {
    public:
        ViewerDataExt(
            igl::opengl::glfw::Viewer* _viewer,
            const Eigen::RowVector3d& color);
        virtual ~ViewerDataExt();
        ViewerData& data() { return m_viewer->data_list[data_id]; }

        virtual void recolor() {}
        virtual inline bool is_mesh() { return false; }
        virtual inline bool is_scalar_field() { return false; }
        virtual inline bool is_vector_field() { return false; }
        virtual inline bool is_com() { return false; }

        /// only used in graph data for now
        virtual void update_vertex_data() {}
        virtual bool visibility() { return data().show_overlay; }
        virtual void visibility(const bool show) { data().show_overlay = show; }

        igl::opengl::glfw::Viewer* m_viewer;
        size_t data_id;

        Eigen::RowVector3d m_color;

        bool show_vertex_data;

        // for scalar field
        bool m_use_log_scale;
        int m_num_isolines;

        // for vector field
        bool m_normalized;
        double m_scaling;

    protected:
        Eigen::MatrixXd set_vertices(const Eigen::MatrixXd& V);

        void set_points(const Eigen::MatrixXd& V, const Eigen::MatrixXd& color);

        void set_edges(
            const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& E,
            const Eigen::MatrixXd& color);

        void set_faces(
            const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& F,
            const Eigen::MatrixXd& color);

        void set_vector_field(
            const Eigen::MatrixXd& V,
            const Eigen::MatrixXd& F,
            const Eigen::MatrixXd& color);
    };

    class MeshData : public ViewerDataExt {
    public:
        MeshData(
            igl::opengl::glfw::Viewer* _viewer,
            const Eigen::RowVector3d& color);

        void set_mesh(
            const Eigen::MatrixXd& V,
            const Eigen::MatrixXi& E,
            const Eigen::MatrixXi& F);

        void set_vertex_data(
            const MatrixXb& data, const Eigen::VectorXi& vertex_type);
        void update_vertex_data() override;
        void update_vertices(const Eigen::MatrixXd& V);

        virtual bool visibility() override { return data().show_overlay; }
        virtual void visibility(const bool show) override
        {
            data().show_overlay = show;
            data().show_lines = show;
        }

        void recolor() override;
        inline bool is_mesh() override { return true; }
        Eigen::MatrixXd mV;
        Eigen::MatrixXi mE;
        Eigen::MatrixXi mF;

        Eigen::VectorXi m_vertex_type;

        Eigen::RowVector3d m_edge_color;
        Eigen::RowVector3d m_static_color;
        Eigen::RowVector3d m_kinematic_color;

        std::vector<std::string> vertex_data_labels;
    };

    class VectorFieldData : public ViewerDataExt {
    public:
        VectorFieldData(
            igl::opengl::glfw::Viewer* _viewer,
            const Eigen::RowVector3d& color);

        virtual inline bool is_vector_field() override { return true; }

        void recolor() override;
        virtual void
        set_vector_field(const Eigen::MatrixXd& V, const Eigen::MatrixXd& F);

        void
        update_vector_field(const Eigen::MatrixXd& V, const Eigen::MatrixXd& F);

        Eigen::MatrixXd mF;
        Eigen::MatrixXd mV;
    };

    class CoMData : public VectorFieldData {
    public:
        CoMData(igl::opengl::glfw::Viewer* _viewer);

        inline bool is_vector_field() override { return false; }
        inline bool is_com() override { return true; }

        void set_coms(const ipc::rigid::PosesD& poses);

        void set_vector_field(
            const Eigen::MatrixXd& V, const Eigen::MatrixXd& F) override;
    };

    class ScalarFieldData : public ViewerDataExt {
    public:
        ScalarFieldData(
            igl::opengl::glfw::Viewer* _viewer,
            const Eigen::RowVector3d& inf_color,
            const Eigen::RowVector3d& bg_color);

        void recolor() override;
        void set_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

        void set_vertex_data(const Eigen::MatrixXd& vertex_data);

        inline bool is_scalar_field() override { return true; }
        bool visibility() override { return is_visible; }
        void visibility(const bool show) override;

    protected:
        void draw_mesh(const Eigen::VectorXd& fx);
        void draw_isolines(const Eigen::VectorXd& fx);
        Eigen::VectorXd inf_to_max(const Eigen::VectorXd& fx);

        Eigen::MatrixXi mF;
        Eigen::MatrixXd mV;

        Eigen::MatrixXd mIsoV;
        Eigen::MatrixXi mIsoE;

        Eigen::MatrixXd m_vertex_data;
        Eigen::MatrixXd m_vertex_clean_data;

        igl::ColorMapType m_colormap_type;
        Eigen::RowVector3d m_base_color;
        bool is_visible;
    };

} // namespace opengl
} // namespace igl
