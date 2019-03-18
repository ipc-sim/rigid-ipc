#include <igl/serialize.h>
#include <opt/minimize.hpp>

namespace ccd {
namespace io {
    void write_opt_results(const std::string filename,
        const ccd::opt::OptimizationResult& results,
        const std::vector<Eigen::MatrixX2d>& u_history,
        const std::vector<double>& f_history,
        const std::vector<double>& g_history);

    void read_opt_results(const std::string fname,
        ccd::opt::OptimizationResult& results,
        std::vector<Eigen::MatrixX2d>& u_history,
        std::vector<double>& f_history, std::vector<double>& g_history);

} // namespace io
} // namespace ccd

namespace igl {
namespace serialization {
    // the `template <>` is essential
    template <>
    inline void serialize(
        const ccd::opt::OptimizationResult& obj, std::vector<char>& buffer)
    {
        ::igl::serialize(obj.x, std::string("x"), buffer);
        ::igl::serialize(obj.fun, std::string("fun"), buffer);
        ::igl::serialize(obj.success, std::string("success"), buffer);
        ::igl::serialize(obj.finished, std::string("finished"), buffer);
        ::igl::serialize(obj.method, std::string("method"), buffer);
    }
    template <>
    inline void deserialize(
        ccd::opt::OptimizationResult& obj, const std::vector<char>& buffer)
    {
        ::igl::deserialize(obj.x, std::string("x"), buffer);
        ::igl::deserialize(obj.fun, std::string("fun"), buffer);
        ::igl::deserialize(obj.success, std::string("success"), buffer);
        ::igl::deserialize(obj.finished, std::string("finished"), buffer);
        ::igl::deserialize(obj.method, std::string("method"), buffer);
    }
} // namespace serialization
} // namespace igl
