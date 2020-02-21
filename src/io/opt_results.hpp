#include <igl/serialize.h>
#include <opt/optimization_results.hpp>

namespace ccd {
namespace io {
    void write_opt_results(
        const std::string filename,
        const ccd::opt::OptimizationResults& results,
        const std::vector<Eigen::MatrixX2d>& u_history,
        const std::vector<double>& f_history,
        const std::vector<double>& g_history);

    void read_opt_results(
        const std::string fname,
        ccd::opt::OptimizationResults& results,
        std::vector<Eigen::MatrixX2d>& u_history,
        std::vector<double>& f_history,
        std::vector<double>& g_history);

} // namespace io
} // namespace ccd

namespace igl {
namespace serialization {
    // the `template <>` is essential
    template <>
    inline void serialize(
        const ccd::opt::OptimizationResults& obj, std::vector<char>& buffer)
    {
        ::igl::serialize(obj.x, std::string("x"), buffer);
        ::igl::serialize(obj.minf, std::string("fun"), buffer);
        ::igl::serialize(obj.success, std::string("success"), buffer);
        ::igl::serialize(obj.finished, std::string("finished"), buffer);
    }
    template <>
    inline void deserialize(
        ccd::opt::OptimizationResults& obj, const std::vector<char>& buffer)
    {
        ::igl::deserialize(obj.x, std::string("x"), buffer);
        ::igl::deserialize(obj.minf, std::string("fun"), buffer);
        ::igl::deserialize(obj.success, std::string("success"), buffer);
        ::igl::deserialize(obj.finished, std::string("finished"), buffer);
    }
} // namespace serialization
} // namespace igl
