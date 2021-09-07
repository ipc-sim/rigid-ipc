#include "serialize_json.hpp"

#include <logger.hpp>
#include <utils/eigen_ext.hpp>

namespace ipc::rigid {

template <>
void from_json<bool>(
    const nlohmann::json& json, Eigen::Matrix<bool, Eigen::Dynamic, 1>& vector)
{
    typedef std::vector<bool> L;
    L list = json.get<L>();
    vector.resize(int(list.size()));
    for (size_t i = 0; i < list.size(); ++i) {
        vector[int(i)] = list[i];
    }
}

nlohmann::json
to_json_string(const Eigen::MatrixXd& matrix, const std::string& format)
{
    size_t num_rows = size_t(matrix.rows());
    size_t num_cols = size_t(matrix.cols());
    std::vector<std::vector<std::string>> vec(num_rows);
    for (size_t i = 0; i < num_rows; ++i) {
        std::vector<std::string> row(num_cols);
        for (size_t j = 0; j < num_cols; ++j) {
            std::string f = "{:" + format + "}";
            row[j] = fmt::format(f, matrix(long(i), long(j)));
        }
        vec[i] = row;
    }
    return nlohmann::json(vec);
}

} // namespace ipc::rigid
