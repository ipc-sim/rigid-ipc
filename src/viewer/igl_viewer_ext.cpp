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
        set_vertices(V);
        set_points(V, color);
        set_edges(V, E, color);

        mE = E;
        mV = V;
        m_color = color;

    }

    void ViewerDataExt::update_graph(const Eigen::MatrixXd& V){
        set_vertices(V);
        set_points(V, m_color);
        set_edges(V, mE, m_color);
    }


    void ViewerDataExt::set_vertices(const Eigen::MatrixXd& V)
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
