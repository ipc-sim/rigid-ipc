#include "opt_results.hpp"

#include <fstream>
#include <iomanip> // std::setw

#include <io/write_scene.hpp>

namespace ccd {
namespace io {

    void write_opt_results(const std::string fname,
        const ccd::opt::OptimizationResult& results,
        const std::vector<Eigen::MatrixX2d>& u_history,
        const std::vector<double>& f_history,
        const std::vector<double>& g_history)
    {
        igl::serialize(results, "results", fname.c_str(), true);
        igl::serialize(u_history, "u_history", fname.c_str());
        igl::serialize(f_history, "f_history", fname.c_str());
        igl::serialize(g_history, "g_history", fname.c_str());
    }

    void read_opt_results(const std::string fname,
        ccd::opt::OptimizationResult& results,
        std::vector<Eigen::MatrixX2d>& u_history,
        std::vector<double>& f_history,
        std::vector<double>& g_history)
    {
        igl::deserialize(results, "results", fname.c_str());
        igl::deserialize(u_history, "u_history", fname.c_str());
        igl::deserialize(f_history, "f_history", fname.c_str());
        igl::deserialize(g_history, "g_history", fname.c_str());
    }

} // namespace io
} // namespace ccd
