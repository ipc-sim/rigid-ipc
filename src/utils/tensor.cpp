#include "tensor.hpp"

namespace ipc::rigid {

Eigen::SparseMatrix<double> multiply(
    const Eigen::VectorXd& x, const std::vector<Eigen::SparseMatrix<double>>& y)
{
    assert(x.rows() == int(y.size()));
    int num_rows, num_cols;
    num_rows = y[0].rows();
    num_cols = y[0].cols();

    Eigen::SparseMatrix<double> m;
    m.resize(num_rows, num_cols);
    m.setZero();

    for (size_t i = 0; i < y.size(); ++i) {
        assert(y[i].rows() == num_rows);
        assert(y[i].cols() == num_cols);
        m += x[int(i)] * y[i];
    }
    return m;
}

} // namespace ipc::rigid
