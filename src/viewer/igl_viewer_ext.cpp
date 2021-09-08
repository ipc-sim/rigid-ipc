#include "igl_viewer_ext.hpp"

#include <igl/colormap.h>
#include <igl/isolines.h>
#include <igl/slice.h>
#include <igl/slice_mask.h>

using namespace ipc;

namespace igl {
namespace opengl {
    ViewerDataExt::ViewerDataExt(
        igl::opengl::glfw::Viewer* _viewer, const Eigen::RowVector3d& color)
        : m_viewer(_viewer)
        , m_color(color)
        , show_vertex_data(false)
        , m_use_log_scale(true)
        , m_normalized(false)
        , m_scaling(1.0)
    {
        data_id = m_viewer->data_list.size() - 1;
    }
    ViewerDataExt::~ViewerDataExt() {}

    Eigen::MatrixXd ViewerDataExt::set_vertices(const Eigen::MatrixXd& V)
    {
        Eigen::MatrixXd V_temp;

        // If P only has two columns, pad with a column of zeros
        if (V.cols() == 2) {
            V_temp = Eigen::MatrixXd::Zero(V.rows(), 3);
            V_temp.block(0, 0, V.rows(), 2) = V;
        } else {
            V_temp = V;
        }
        data().set_vertices(V_temp);

        Eigen::MatrixXd nz = Eigen::MatrixXd::Zero(V.rows(), 3);
        nz.col(2).setConstant(1.0);
        data().set_normals(nz);
        return V_temp;
    }

    void ViewerDataExt::set_points(
        const Eigen::MatrixXd& V, const Eigen::MatrixXd& color)
    {
        Eigen::MatrixXd V_temp;

        // If P only has two columns, pad with a column of zeros
        if (V.cols() == 2) {
            V_temp = Eigen::MatrixXd::Zero(V.rows(), 3);
            V_temp.block(0, 0, V.rows(), 2) = V;
        } else {
            V_temp = V;
        }
        data().set_points(V_temp, color);
    }
    void ViewerDataExt::set_edges(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXd& color)
    {
        Eigen::MatrixXd V_temp;

        // If P only has two columns, pad with a column of zeros
        if (V.cols() == 2) {
            V_temp = Eigen::MatrixXd::Zero(V.rows(), 3);
            V_temp.block(0, 0, V.rows(), 2) = V;
        } else {
            V_temp = V;
        }
        data().set_edges(V_temp, E, color);
    }
    void ViewerDataExt::set_faces(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::MatrixXd& color)
    {
        data().V = Eigen::MatrixXd(0, 3);
        data().F = Eigen::MatrixXi(0, 3);
        data().set_mesh(V, F);
        data().set_colors(color);
        data().show_lines = false;
    }

    ///////////////////////////////////////////////////////////////////////////
    /// \brief MeshData::MeshData
    /// \param _viewer
    ////////////////////////////////////////////////////////////////////////
    MeshData::MeshData(
        igl::opengl::glfw::Viewer* _viewer, const Eigen::RowVector3d& color)
        : ViewerDataExt(_viewer, color)
    {
        m_edge_color = Eigen::RowVector3d::Zero();                    // #000000
        m_static_color = Eigen::RowVector3d(0xB3, 0xB3, 0xB3) / 0xFF; // #B3B3B3
        m_kinematic_color =
            Eigen::RowVector3d(0xFF, 0x80, 0x00) / 0xFF; // #FF8000
    }

    void MeshData::set_mesh(
        const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::MatrixXi& F)
    {
        data().clear();
        mV = set_vertices(V);
        set_points(V, F.size() ? m_edge_color : m_color);
        set_edges(V, E, F.size() ? m_edge_color : m_color);
        if (F.size()) {
            set_faces(V, F, m_color);
        }

        mE = E;
        mF = F;
    }

    void MeshData::update_vertices(const Eigen::MatrixXd& V)
    {
        mV = set_vertices(V);
        recolor();
        data().labels_positions = mV;
    }

    void MeshData::recolor()
    {
        assert(mV.rows() == m_vertex_type.rows());

        Eigen::MatrixXd vertex_colors(mV.rows(), 3);
        Eigen::MatrixXd edge_vertex_colors(mE.rows(), 3);
        Eigen::MatrixXd face_vertex_colors(mV.rows(), 3);

        auto type_color = [&](int type) {
            return type == 0 ? m_static_color
                             : (type == 1 ? m_kinematic_color : m_color);
        };

        for (int i = 0; i < mV.rows(); i++) {
            vertex_colors.row(i) =
                mF.size() ? m_edge_color : type_color(m_vertex_type(i));
            face_vertex_colors.row(i) = type_color(m_vertex_type(i));
        }

        for (int i = 0; i < mE.rows(); i++) {
            edge_vertex_colors.row(i) = mF.size()
                ? m_edge_color
                : type_color(std::min(
                      m_vertex_type(mE(i, 0)), m_vertex_type(mE(i, 1))));
        }

        set_points(mV, vertex_colors);
        if (mE.size()) {
            set_edges(mV, mE, edge_vertex_colors);
        }
        if (mF.size()) {
            set_faces(mV, mF, face_vertex_colors);
        }
    }

    void MeshData::set_vertex_data(
        const MatrixXb& vtx_data, const Eigen::VectorXi& vertex_type)
    {
        assert(vtx_data.rows() == mV.rows());

        vertex_data_labels.clear();
        vertex_data_labels.resize(size_t(vtx_data.rows()));
        for (size_t i = 0; i < vtx_data.rows(); i++) {
            std::string data_i = std::to_string(i) + ":";
            for (size_t j = 0; j < vtx_data.cols(); j++) {
                if (vtx_data(i, j)) {
                    data_i += std::to_string(j) + ",";
                }
            }
            vertex_data_labels[i] = data_i.substr(0, data_i.size() - 1);
        }
        show_vertex_data = true;
        data().labels_positions = mV;
        data().labels_strings = vertex_data_labels;

        m_vertex_type = vertex_type;
        recolor();
    }

    void MeshData::update_vertex_data()
    {
        if (show_vertex_data) {
            data().labels_positions = mV;
            data().labels_strings = vertex_data_labels;
        } else {
            data().labels_positions.resize(0, 0);
        }
    }

    ////////////////////////////////////////////////////////////////////////
    /// \brief VectorFieldData::VectorFieldData
    /// \param _viewer
    ////////////////////////////////////////////////////////////////////////
    VectorFieldData::VectorFieldData(
        igl::opengl::glfw::Viewer* _viewer, const Eigen::RowVector3d& color)
        : ViewerDataExt(_viewer, color)
    {
    }
    void VectorFieldData::set_vector_field(
        const Eigen::MatrixXd& V, const Eigen::MatrixXd& F)
    {
        Eigen::MatrixXd f = F;
        if (f.size() == 0) {
            f.resizeLike(V);
            f.setZero();
        }
        if (m_normalized) {
            f.rowwise().normalize();
        }
        f *= m_scaling;

        data().lines.resize(0, 9);
        data().add_edges(V, V + f, m_color);
        mV = V;
        mF = F;
    }

    void VectorFieldData::update_vector_field(
        const Eigen::MatrixXd& V, const Eigen::MatrixXd& F)
    {
        set_vector_field(V, F);
    }

    void VectorFieldData::recolor()
    {
        if (mF.size() > 0) {
            set_vector_field(mV, mF);
        }
    }

    ////////////////////////////////////////////////////////////////////////
    /// \brief CoMData::CoMData
    /// \param _viewer
    ////////////////////////////////////////////////////////////////////////
    CoMData::CoMData(igl::opengl::glfw::Viewer* _viewer)
        : VectorFieldData(_viewer, Eigen::RowVector3d::Zero())
    {
    }

    void CoMData::set_coms(const ipc::rigid::PosesD& poses)
    {
        int dim = poses.size() ? poses[0].dim() : 0;
        Eigen::MatrixXd com(dim * poses.size(), dim);
        Eigen::MatrixXd principle_axes(dim * poses.size(), dim);

        for (int i = 0; i < poses.size(); i++) {
            MatrixMax3d R = poses[i].construct_rotation_matrix();
            for (int j = 0; j < dim; j++) {
                com.row(dim * i + j) = poses[i].position;
                VectorMax3d axis = VectorMax3d::Zero(dim);
                axis(j) = 1;
                principle_axes.row(dim * i + j) = R * axis;
            }
        }

        this->set_vector_field(com, principle_axes);
    }

    void CoMData::set_vector_field(
        const Eigen::MatrixXd& V, const Eigen::MatrixXd& F)
    {
        Eigen::MatrixXd f = F;
        if (f.size() == 0) {
            f.resizeLike(V);
            f.setZero();
        }
        if (m_normalized) {
            f.rowwise().normalize();
        }
        f *= m_scaling;

        data().lines.resize(0, 9);

        int dim = V.cols();
        Eigen::MatrixXd colors(V.rows(), 3);
        for (int i = 0; i < colors.rows(); i += dim) {
            colors.row(i + 0) =
                Eigen::RowVector3d(0xFF, 0x26, 0x00) / 0xFF; // #FF2600
            colors.row(i + 1) =
                Eigen::RowVector3d(0x00, 0xF9, 0x00) / 0xFF; // #00F900
            if (dim == 3) {
                colors.row(i + 2) =
                    Eigen::RowVector3d(0x03, 0x29, 0xD0) / 0xFF; // #0329D0
            }
        }

        data().add_edges(V, V + f, colors);
        mV = V;
        mF = F;
    }

    ////////////////////////////////////////////////////////////////////////
    /// \brief ScalarFieldData::ScalarFieldData
    /// \param _viewer
    ////////////////////////////////////////////////////////////////////////
    ScalarFieldData::ScalarFieldData(
        igl::opengl::glfw::Viewer* _viewer,
        const Eigen::RowVector3d& inf_color,
        const Eigen::RowVector3d& bg_color)
        : ViewerDataExt(_viewer, inf_color)
        , m_colormap_type(igl::COLOR_MAP_TYPE_JET)
        , m_base_color(bg_color)
        , is_visible(false)
    {
        m_num_isolines = 10;
    }

    void ScalarFieldData::set_mesh(
        const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
    {
        mV = V;
        mF = F;
        data().clear();
        data().set_mesh(V, F);
        data().set_colors(m_base_color);
        data().show_lines = false;
        data().show_faces = true;
        data().show_overlay = true;
    }

    void ScalarFieldData::set_vertex_data(const Eigen::MatrixXd& fx)
    {
        m_vertex_data = fx;
        m_vertex_clean_data = inf_to_max(fx);
        recolor();
    }

    void ScalarFieldData::recolor()
    {
        if (m_vertex_data.size() == 0 || !is_visible) {
            return;
        }
        Eigen::VectorXd fx = m_vertex_clean_data;
        if (m_use_log_scale) {
            fx = fx.array().log();
        }

        draw_mesh(fx);
        draw_isolines(fx);
    }

    void ScalarFieldData::visibility(const bool show)
    {
        static bool _show_faces = data().show_faces;
        static bool _show_overlay = data().show_overlay;

        is_visible = show;
        if (is_visible) {
            data().show_faces = _show_faces;
            data().show_overlay = _show_overlay;
        } else {
            _show_faces = data().show_faces;
            _show_overlay = data().show_overlay;
            data().show_faces = false;
            data().show_overlay = false;
        }
    }

    void ScalarFieldData::draw_isolines(const Eigen::VectorXd& fx)
    {
        igl::isolines(mV, mF, fx, m_num_isolines, mIsoV, mIsoE);
        set_edges(mIsoV, mIsoE, m_color);
    }

    void ScalarFieldData::draw_mesh(const Eigen::VectorXd& fx)
    {
        Eigen::MatrixXd c;
        igl::colormap(igl::COLOR_MAP_TYPE_JET, fx, /*normalize=*/true, c);

        for (int i = 0; i < m_vertex_data.rows(); ++i) {
            if (std::isinf(m_vertex_data(i))) {
                c.row(i) = m_color;
            }
        }
        data().set_colors(c);
    }
    Eigen::VectorXd ScalarFieldData::inf_to_max(const Eigen::VectorXd& fx)
    {
        auto inf_to_neg = [](const double x) {
            if (std::isinf(x)) {
                return -1.0;
            }
            return x;
        };

        Eigen::VectorXd fx_clean = fx.unaryExpr(inf_to_neg);
        double max_coeff = fx_clean.maxCoeff();
        // neg to max_coeff
        fx_clean = (fx_clean.array() < 0).select(max_coeff, fx_clean);
        return fx_clean;
    }
} // namespace opengl
} // namespace igl
