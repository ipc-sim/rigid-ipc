#include <state.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

int main(int argc, char* argv[])
{
    std::string out_dir = std::string(DATA_OUTPUT_DIR);
    if (argc >= 2) {
        out_dir = argv[1];
    }

    ccd::State state;
    state.solver_settings.max_iter = 1000;
    state.solver_settings.method = ccd::opt::NCP;
    state.solver_settings.lcp_solver = ccd::opt::LCP_MOSEK;

    Eigen::MatrixX2d vertices;
    Eigen::MatrixX2i edges;
    Eigen::MatrixX2d displacements;


    state.load_scene(vertices, edges, displacements);
    state.optimize_displacements();

    return 0;
}
