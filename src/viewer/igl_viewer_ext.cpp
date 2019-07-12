#include "igl_viewer_ext.hpp"
namespace igl {
namespace opengl {
    ViewerDataExt::ViewerDataExt(igl::opengl::glfw::Viewer* _viewer)
        : m_viewer(_viewer)
        , show_vertex_data(false)
    {
        data_id = m_viewer->data_list.size() - 1;
    }

    void ViewerDataExt::set_graph(const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E,
        const Eigen::RowVector3d& color)
    {
        mV = set_vertices(V);
        set_points(V, color);
        set_edges(V, E, color);

        mE = E;
        m_color = color;
    }

    void ViewerDataExt::update_graph(const Eigen::MatrixXd& V)
    {
        mV = set_vertices(V);
        set_points(V, m_color);
        set_edges(V, mE, m_color);
        data().labels_positions = mV;
    }

    void ViewerDataExt::set_vector_field(const Eigen::MatrixXd& V,
        const Eigen::MatrixXd& F,
        const Eigen::RowVector3d& color)
    {
        data().lines.resize(0, 9);
        data().add_edges(V, V + F, color);
        mV = V;
        mF = F;
        m_color = color;
    }

    void ViewerDataExt::update_vector_field(
        const Eigen::MatrixXd& V, const Eigen::MatrixXd& F)
    {
        set_vector_field(V, F, m_color);
    }
    void ViewerDataExt::recolor(){
        if (mF.size() > 0){
            set_vector_field(mV, mF, m_color);
        }
        if (mE.size() > 0){
            set_points(mV, m_color);
            set_edges(mV, mE, m_color);
        }
    }

    void ViewerDataExt::set_vertex_data(
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

    void ViewerDataExt::update_vertex_data(){
        if (show_vertex_data){
            data().labels_positions = mV;
            data().labels_strings = vertex_data_labels;
        }
        else {
            data().labels_positions.resize(0,0);
        }
    }

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

} // namespace opengl
} // namespace igl
