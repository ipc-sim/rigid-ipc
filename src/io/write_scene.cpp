#include "io/write_scene.hpp"

#include <fstream>
#include <iomanip> // std::setw

namespace ccd {
namespace io {
    template <typename type>
    std::vector<std::vector<type>> to_vec2(
        const Eigen::Matrix<type, Eigen::Dynamic, Eigen::Dynamic>& matrix)
    {
        size_t num_rows = size_t(matrix.rows());
        size_t num_cols = size_t(matrix.cols());
        std::vector<std::vector<type>> vec(num_rows);
        for (size_t i = 0; i < num_rows; ++i) {
            std::vector<type> row(num_cols);
            for (size_t j = 0; j < num_cols; ++j) {
                row[j] = matrix(long(i), long(j));
            }
            vec[i] = row;
        }
        return vec;
    }

    void write_scene(const std::string filename,
        const Eigen::MatrixXd& vertices, const Eigen::MatrixX2i& edges,
        const Eigen::MatrixXd& displacements)
    {
        using nlohmann::json;
        json scene = write_scene(vertices, edges, displacements);

        std::ofstream o(filename);
        o << std::setw(4) << scene << std::endl;
    }

    nlohmann::json write_scene(const Eigen::MatrixXd& vertices,
        const Eigen::MatrixX2i& edges, const Eigen::MatrixXd& displacements)
    {
        using nlohmann::json;
        json scene;
        scene["vertices"] = to_vec2<double>(vertices);
        scene["displacements"] = to_vec2<double>(displacements);
        scene["edges"] = to_vec2<int>(edges);
        return scene;
    }
} // namespace io
} // namespace ccd
