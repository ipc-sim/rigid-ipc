#include "igl_viewer_ext.hpp"
namespace igl {
namespace opengl {
    ViewerDataExt::ViewerDataExt(ViewerData& data)
        : m_data(data)
    {
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

    void ViewerDataExt::set_vertex_data(const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> data){
        assert(data.rows() == mV.rows());

        vertex_data_labels.clear();
        vertex_data_labels.resize(size_t(data.rows()));
        for (uint i=0; i < data.rows(); i++){
            std::string data_i = std::to_string(i) +":";
            for (uint j=0; j < data.cols(); j++){
                if (data(i,j)){
                    data_i += std::to_string(j) + ",";
                }
            }
            vertex_data_labels[i] = data_i.substr(0,data_i.size()-1);
        }

        m_data.labels_positions = mV;
        m_data.labels_strings = vertex_data_labels;

    }

    void ViewerDataExt::update_graph(const Eigen::MatrixXd& V){
        mV = set_vertices(V);
        set_points(V, m_color);
        set_edges(V, mE, m_color);
        m_data.labels_positions = mV;
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
        m_data.set_vertices(V_temp);

        Eigen::MatrixXd nz = Eigen::MatrixXd::Zero(V.rows(), 3);
        nz.col(2).setConstant(1.0);
        m_data.set_normals(nz);
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
        m_data.set_points(V_temp, color);
    }
    void ViewerDataExt::set_edges(const Eigen::MatrixXd& V,
        const Eigen::MatrixXi& E, const Eigen::MatrixXd& color)
    {
        Eigen::MatrixXd V_temp;

        // If P only has two columns, pad with a column of zeros
        if (V.cols() == 2) {
            V_temp = Eigen::MatrixXd::Zero(V.rows(), 3);
            V_temp.block(0, 0, V.rows(), 2) = V;
        } else {
            V_temp = V;
        }
        m_data.set_edges(V_temp, E, color);
    }

} // namespace opengl
} // namespace igl
