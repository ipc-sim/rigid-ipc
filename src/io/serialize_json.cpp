#include "serialize_json.hpp"
#include <utils/eigen_ext.hpp>

namespace ccd {
namespace io {

template <>
void from_json<bool>(const nlohmann::json& json,
    Eigen::Matrix<bool, Eigen::Dynamic, 1>& vector)
{
    typedef std::vector<bool> L;
    L list = json.get<L>();
    vector.resize(int(list.size()));
    for (size_t i = 0; i < list.size(); ++i) {
        vector[int(i)] = list[i];
    }
}

} // namespace io
} // namespace ccd
