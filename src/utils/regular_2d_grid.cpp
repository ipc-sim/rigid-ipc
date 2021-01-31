#include "regular_2d_grid.hpp"

#include <array>
#include <vector>

namespace ipc::rigid {
void regular_2d_grid(
    const int n, const bool tri, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{

    V.resize(n * n, 2);
    F.resize((n - 1) * (n - 1) * 2, 3);
    const double delta = 1. / (n - 1.);
    std::vector<int> map(size_t(n * n), -1);

    int index = 0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (tri && i + j >= n) {
                continue;
            }
            map[size_t(i + j * n)] = index;
            V.row(index) << i * delta, j * delta;
            ++index;
        }
    }

    V.conservativeResize(index, 2);

    std::array<int, 3> tmp;

    index = 0;
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - 1; ++j) {
            tmp = { { map[i + j * n], map[i + 1 + j * n],
                      map[i + (j + 1) * n] } };
            if (tmp[0] >= 0 && tmp[1] >= 0 && tmp[2] >= 0) {
                F.row(index) << tmp[0], tmp[1], tmp[2];
                ++index;
            }

            tmp = { { map[i + 1 + j * n], map[i + 1 + (j + 1) * n],
                      map[i + (j + 1) * n] } };
            if (tmp[0] >= 0 && tmp[1] >= 0 && tmp[2] >= 0) {
                F.row(index) << tmp[0], tmp[1], tmp[2];
                ++index;
            }
        }
    }

    F.conservativeResize(index, 3);
}

void regular_2d_grid(
    const int num_cols,
    const int num_rows,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F)
{

    V.resize(num_cols * num_rows, 2);
    F.resize((num_cols - 1) * (num_rows - 1) * 2, 3);
    const double delta_row = 1. / (num_rows - 1.);
    const double delta_col = 1. / (num_cols - 1.);
    std::vector<int> map(size_t(num_cols * num_rows), -1);

    int index = 0;
    for (int i = 0; i < num_cols; ++i) {
        for (int j = 0; j < num_rows; ++j) {
            map[size_t(i + j * num_cols)] = index;
            V.row(index) << i * delta_col, j * delta_row;
            ++index;
        }
    }

    V.conservativeResize(index, 2);

    std::array<int, 3> tmp;

    index = 0;
    for (int i = 0; i < num_cols - 1; ++i) {
        for (int j = 0; j < num_rows - 1; ++j) {
            tmp = { { map[i + j * num_cols], map[i + 1 + j * num_cols],
                      map[i + (j + 1) * num_cols] } };
            if (tmp[0] >= 0 && tmp[1] >= 0 && tmp[2] >= 0) {
                F.row(index) << tmp[0], tmp[1], tmp[2];
                ++index;
            }

            tmp = { { map[i + 1 + j * num_cols],
                      map[i + 1 + (j + 1) * num_cols],
                      map[i + (j + 1) * num_cols] } };
            if (tmp[0] >= 0 && tmp[1] >= 0 && tmp[2] >= 0) {
                F.row(index) << tmp[0], tmp[1], tmp[2];
                ++index;
            }
        }
    }

    F.conservativeResize(index, 3);
}
} // namespace ipc::rigid
