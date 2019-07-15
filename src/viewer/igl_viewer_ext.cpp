#include "igl_viewer_ext.hpp"

#include <igl/colormap.h>
#include <igl/slice.h>
#include <igl/slice_mask.h>

#include <utils/eigen_ext.hpp>

namespace igl {
namespace opengl {
    ViewerDataExt::ViewerDataExt(igl::opengl::glfw::Viewer* _viewer)
        : m_viewer(_viewer)
        , show_vertex_data(false)
        , m_log_scale(true)
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
    void ViewerDataExt::set_edges(const Eigen::MatrixXd& V,
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

    ///////////////////////////////////////////////////////////////////////////
    /// \brief GraphData::GraphData
    /// \param _viewer
    ////////////////////////////////////////////////////////////////////////
    GraphData::GraphData(igl::opengl::glfw::Viewer* _viewer)
        : ViewerDataExt(_viewer)

    {
    }

    void GraphData::set_graph(const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::RowVector3d& color)
    {
        mV = set_vertices(V);
        set_points(V, color);
        set_edges(V, E, color);

        mE = E;
        m_color = color;
    }

    void GraphData::update_graph(const Eigen::MatrixXd& V)
    {
        mV = set_vertices(V);
        set_points(V, m_color);
        set_edges(V, mE, m_color);
        data().labels_positions = mV;
    }

    void GraphData::recolor()
    {
        if (mE.size() > 0) {
            set_points(mV, m_color);
            set_edges(mV, mE, m_color);
        }
    }
    void GraphData::set_vertex_data(
        const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> vtx_data)
    {
        assert(vtx_data.rows() == mV.rows());

        vertex_data_labels.clear();
        vertex_data_labels.resize(size_t(vtx_data.rows()));
        for (uint i = 0; i < vtx_data.rows(); i++) {
            std::string data_i = std::to_string(i) + ":";
            for (uint j = 0; j < vtx_data.cols(); j++) {
                if (vtx_data(i, j)) {
                    data_i += std::to_string(j) + ",";
                }
            }
            vertex_data_labels[i] = data_i.substr(0, data_i.size() - 1);
        }
        show_vertex_data = true;
        data().labels_positions = mV;
        data().labels_strings = vertex_data_labels;
    }

    void GraphData::update_vertex_data()
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
    VectorFieldData::VectorFieldData(igl::opengl::glfw::Viewer* _viewer)
        : ViewerDataExt(_viewer)
    {
    }
    void VectorFieldData::set_vector_field(const Eigen::MatrixXd& V,
        const Eigen::MatrixXd& F,
        const Eigen::RowVector3d& color)
    {
        data().lines.resize(0, 9);
        data().add_edges(V, V + F, color);
        mV = V;
        mF = F;
        m_color = color;
    }

    void VectorFieldData::update_vector_field(
        const Eigen::MatrixXd& V, const Eigen::MatrixXd& F)
    {
        set_vector_field(V, F, m_color);
    }

    void VectorFieldData::recolor()
    {
        if (mF.size() > 0) {
            set_vector_field(mV, mF, m_color);
        }
    }

    ////////////////////////////////////////////////////////////////////////
    /// \brief ScalarFieldData::ScalarFieldData
    /// \param _viewer
    ////////////////////////////////////////////////////////////////////////
    ScalarFieldData::ScalarFieldData(igl::opengl::glfw::Viewer* _viewer)
        : ViewerDataExt(_viewer)
        , m_colormap_type(igl::COLOR_MAP_TYPE_JET)

    {
    }

    void ScalarFieldData::set_mesh(const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& F,
        const Eigen::RowVector3d& base_color,
        const Eigen::RowVector3d& inf_color)
    {
        mV = V;
        mF = F;
        data().set_mesh(V, F);
        data().set_colors(base_color);
        data().show_lines = false;
        m_base_color = base_color;
        m_color = inf_color;
    }

    void ScalarFieldData::set_vertex_data(const Eigen::MatrixXd& fx)
    {
        m_vertex_data = fx;
        recolor();
    }

    void ScalarFieldData::recolor()
    {
        Eigen::VectorXb is_finite;
        is_finite = m_vertex_data
                        .unaryExpr([](const double x) {
                            if (std::isinf(x)) {
                                return 0.0;
                            }
                            return 1.0;
                        })
                        .cast<bool>();

        Eigen::VectorXd fx_clean
            = igl::slice_mask(m_vertex_data, is_finite.array(), 1);
        if (m_log_scale) {
            fx_clean = fx_clean.array().log();
        }

        Eigen::MatrixXd c;
        igl::colormap(igl::COLOR_MAP_TYPE_JET, fx_clean, /*normalize=*/true, c);

        Eigen::MatrixXd C(m_vertex_data.rows(), 3);
        int j = 0;
        for (int i = 0; i < m_vertex_data.rows(); ++i) {
            if (std::isinf(m_vertex_data(i))) {
                C.row(i) = m_color;
            } else {
                C.row(i) = c.row(j);
                j += 1;
            }
        }

        data().set_colors(C);
    }
} // namespace opengl
} // namespace igl
