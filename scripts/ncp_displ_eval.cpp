#include <state.hpp>

#include <cmath>
#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

const int NUM_EDGES = 2;
const int NUM_VERTICES = 4;
const int ANGLE_STEPS = 10;

int main(int argc, char* argv[])
{
    using namespace ccd::opt;

    std::string out_dir = std::string(DATA_OUTPUT_DIR);
    if (argc >= 2) {
        out_dir = argv[1];
    }
    spdlog::set_level(spdlog::level::info);

    ccd::State state;
    state.solver_settings.max_iter = 100;
    state.solver_settings.method = ccd::opt::NCP;
    state.volume_epsilon = 1E-6;

    Eigen::MatrixX2d vertices(NUM_VERTICES, 2);
    Eigen::MatrixX2i edges(NUM_EDGES, 2);
    Eigen::MatrixX2d displacements(NUM_VERTICES, 2);

    edges.row(0) << 0, 1;
    edges.row(1) << 2, 3;

    // horizontal fixed edge of length 1.0
    vertices.row(0) << -0.5, 0.0;
    vertices.row(1) << 0.5, 0.0;

    displacements.row(0) << 0.0, 0.0;
    displacements.row(1) << 0.0, 0.0;

    LCPSolver lcp_solvers[] = { LCP_GAUSS_SEIDEL, LCP_MOSEK };

    NcpUpdate update_types[] = { NcpUpdate::G_GRADIENT, NcpUpdate::LINEARIZED };

    double alpha = 0.5;

    Eigen::MatrixX2d expected(NUM_VERTICES, 2);

    std::string summary = out_dir + "/displ_ncp/summary.csv";
    std::ofstream o(summary);
    o << "theta, dy, f(x*), f(e*), sum g(x*), sum g(e*), num_it, x_feasible, "
         "x_optimal, e_feasible, path"
      << std::endl;

    for (int i = 0; i < ANGLE_STEPS + 1; i++) {
        double theta = M_PI * i / 2.0 / ANGLE_STEPS;

        vertices.row(2) << -0.5 + alpha, 0.5;
        vertices.row(3) << -0.5 + alpha + cos(theta), sin(theta) + 0.6;

        for (int j = 0; j < 10; j++) {
            double dy = (0.5 + 1E-4) * (10.0 - j) / 10.0 + 1.0 * j / 10.0;

            // vertical fall
            displacements.row(2) << 0.0, -dy;
            displacements.row(3) << 0.0, -dy;

            double edy = (dy - 0.5) / 2.0;
            expected.row(0) << 0.0, -edy;
            expected.row(1) << 0.0, -edy;

            expected.row(2) << 0.0, -dy + edy;
            expected.row(3) << 0.0, -dy;

            Eigen::MatrixXd e = expected;
            e.resize(e.size(), 1);

            state.load_scene(vertices, edges, displacements);
            std::string scene_name = fmt::format(
                "/displ_ncp/test_alpha_{}_theta_{}_dy_{}", alpha, i, j);
            state.save_scene(out_dir + scene_name + ".json");

            for (const auto& update_type : update_types) {
                for (const auto& lcp_solver : lcp_solvers) {
                    state.solver_settings.ncp_update_method = update_type;
                    state.solver_settings.lcp_solver = lcp_solver;
                    std::string filename = fmt::format("{}_{}_{}", scene_name,
                        LCPSolverNames[lcp_solver],
                        NcpUpdateNames[static_cast<int>(update_type)]);

//                    spdlog::info("case {}", filename);
                    state.optimize_displacements(out_dir + filename + ".csv");

                    double fe = state.opt_problem.f(e);
                    double ge_sum = state.opt_problem.g(e).sum();
                    Eigen::MatrixXd u_ = state.opt_results.x;
                    u_.resize(u_.size(), 1);

                    double fx = state.opt_problem.f(u_);
                    double gx_sum = state.opt_problem.g(u_).sum();

                    o << theta << ",";
                    o << dy << ",";
                    o << fx << ",";
                    o << fe << ",";
                    o << gx_sum << ",";
                    o << ge_sum << ",";
                    o << state.u_history.size() << ",";
                    o << int(gx_sum >= 0) << ","; /// x_feasible
                    o << int(fx <= fe) << ",";    /// x_optimal
                    o << int(ge_sum >= 0) << ","; /// e_feasible
                    o << filename << ",";
                    o << std::endl;


                }
            }
        }
        o.flush();
        spdlog::info("outer loop {}", i);
    }
    o.close();

    return 0;
}
