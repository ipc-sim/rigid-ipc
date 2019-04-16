#include <opt/ncp_solver.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <autodiff/autodiff_types.hpp>
#include <autodiff/finitediff.hpp>

// ---------------------------------------------------
// format for logs
// ---------------------------------------------------
std::string fmt_red(const std::string& s)
{
    return "\033[1;31m" + s + "\033[0m";
}
std::string fmt_green(const std::string& s)
{
    return "\033[1;32m" + s + "\033[0m";
}

std::string fmt_bool(const bool s)
{
    if (!s) {
        return fmt_red("F");
    }
    return fmt_green("T");
}
static const Eigen::IOFormat CommaFmt(
    Eigen::StreamPrecision, Eigen::DontAlignCols, ", ", ", ", "", "", "", "");

std::string fmt_eigen(const Eigen::MatrixXd& x, const int precision = 10)
{
    std::stringstream ssx;
    ssx << std::setprecision(precision) << x.format(CommaFmt);
    return ssx.str();
}

// ---------------------------------------------------
// SETUP
// ---------------------------------------------------
static const int NUM_VARS = 2;
static const int NUM_CONSTRAINTS = 2;

// differentiable helpers
template <typename T> using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;

typedef Eigen::Matrix<double, 2, 1> Vector2d;
typedef DScalar1<double, Eigen::Matrix<double, NUM_VARS, 1>> DScalar;
typedef VectorXT<DScalar> DVector;

// ---------------------------------------------------
// CASES
// ---------------------------------------------------

enum class GxCases { LINEAR, QUADRATIC, ABS_VALUE, CIRCLE, MAX, SMOOTH_MAX };
static const char* GxCasesNames[]
    = { "LINEAR", "QUADRATIC", "ABS_VALUE", "CIRCLE", "MAX", "SMOOTH_MAX" };
const GxCases GxCasesAll[] = { GxCases::LINEAR, GxCases::QUADRATIC,
    GxCases::ABS_VALUE, GxCases::CIRCLE, GxCases::MAX, GxCases::SMOOTH_MAX };

template <typename dscalar>
void get_linear_case(Eigen::VectorXd& expected,
    std::function<VectorXT<dscalar>(const Eigen::VectorXd& x)>& g_diff)
{

    using namespace ccd::opt;
    g_diff = [](const Eigen::VectorXd& x) -> VectorXT<dscalar> {
        VectorXT<dscalar> gx(NUM_CONSTRAINTS);
        dscalar x0(0, x[0]);
        dscalar x1(1, x[1]);

        gx(0) = x0;
        gx(1) = x1;
        return gx;
    };
    expected << 0.0, 0.0;
}

template <typename dscalar>
void get_quadratic_case(Eigen::VectorXd& expected,
    std::function<VectorXT<dscalar>(const Eigen::VectorXd& x)>& g_diff)
{

    using namespace ccd::opt;
    g_diff = [](const Eigen::VectorXd& x) -> VectorXT<dscalar> {
        VectorXT<dscalar> gx(NUM_CONSTRAINTS);
        dscalar x0(0, x[0]);
        dscalar x1(1, x[1]);

        gx(0) = 0.04 - x0 * x0;
        gx(1) = 0.09 - x1 * x1;
        return gx;
    };
    expected << -0.2, -0.3;
}

template <typename dscalar>
void get_abs_value_case(Eigen::VectorXd& expected,
    std::function<VectorXT<dscalar>(const Eigen::VectorXd& x)>& g_diff)
{

    using namespace ccd::opt;
    g_diff = [](const Eigen::VectorXd& x) -> VectorXT<dscalar> {
        VectorXT<dscalar> gx(NUM_CONSTRAINTS);
        dscalar x0(0, x[0]);
        dscalar x1(1, x[1]);

        gx(0) = 0.2 - (x0 > 0 ? x0 : -x0);
        gx(1) = 0.3 - (x1 > 0 ? x1 : -x1);
        return gx;
    };
    expected << -0.2, -0.3;
}

template <typename dscalar>
void get_circle_case(Eigen::VectorXd& expected,
    std::function<VectorXT<dscalar>(const Eigen::VectorXd& x)>& g_diff)
{

    using namespace ccd::opt;
    g_diff = [](const Eigen::VectorXd& x) -> VectorXT<dscalar> {
        VectorXT<dscalar> gx(NUM_CONSTRAINTS);
        dscalar x0(0, x[0]);
        dscalar x1(1, x[1]);

        gx(0) = 1.0 - (x0 - 1.0) * (x0 - 1.0);
        gx(1) = 1.0 - (x1 - 2.5) * (x1 - 2.5);
        return gx;
    };
    expected << 0.0, 1.5;
}

template <typename dscalar>
void get_max_case(Eigen::VectorXd& expected,
    std::function<VectorXT<dscalar>(const Eigen::VectorXd& x)>& g_diff)
{

    using namespace ccd::opt;
    g_diff = [](const Eigen::VectorXd& x) -> VectorXT<dscalar> {
        VectorXT<dscalar> gx(NUM_CONSTRAINTS);
        dscalar x0(0, x[0]);
        dscalar x1(1, x[1]);

        auto aux1 = x1 * x1;
        auto aux0 = x0 * x0;
        gx(0) = 1.0 - (aux0 > 1 ? aux0 : dscalar(1.0));
        gx(1) = 1.0 - (aux1 > 1 ? aux1 : dscalar(1.0));
        return gx;
    };
    expected << -1.0, -1.0;
}

template <typename dscalar>
void get_smooth_max_case(Eigen::VectorXd& expected,
    std::function<VectorXT<dscalar>(const Eigen::VectorXd& x)>& g_diff)
{

    using namespace ccd::opt;
    // s.t  (x0 - 1)^2 <= 1, (x1 - 2.5)^2 <= 1
    g_diff = [](const Eigen::VectorXd& x) -> VectorXT<dscalar> {
        VectorXT<dscalar> gx(NUM_CONSTRAINTS);
        dscalar x0(0, x[0]);
        dscalar x1(1, x[1]);

        auto abs_x0 = x0 > 0 ? x0 : -x0;
        auto abs_x1 = x1 > 0 ? x1 : -x1;
        auto min_x0 = abs_x0 <= 1 ? abs_x0 : dscalar(1.0);
        auto min_x1 = abs_x1 <= 1 ? abs_x1 : dscalar(1.0);
        auto aux_x0 = (x0 - x0 / abs_x0 * min_x0);
        auto aux_x1 = (x1 - x1 / abs_x1 * min_x1);
        gx(0) = -aux_x0 * aux_x0;
        gx(1) = -aux_x1 * aux_x1;
        return gx;
    };
    expected << -1.0, -1.0;
}

// ---------------------------------------------------
// Helper Implementation
// ---------------------------------------------------
void eval_ncp(const GxCases gx_case, const ccd::opt::NcpUpdate update_type,
    const ccd::opt::LCPSolver lcp_solver, const bool check_convergence,
    const std::string outfile)
{

    using namespace ccd::opt;
    DiffScalarBase::setVariableCount(size_t(NUM_VARS));

    Eigen::SparseMatrix<double> A(NUM_VARS, NUM_VARS);
    Eigen::VectorXd b(NUM_VARS), expected(NUM_VARS);

    // ------------------------------------------------------------------------
    // PROBLEM SETUP
    // ------------------------------------------------------------------------
    A.setIdentity();
    b << -1, -2.5;

    std::function<DVector(const Eigen::VectorXd& x)> g_diff;

    switch (gx_case) {
    case GxCases::LINEAR:
        get_linear_case(expected, g_diff);
        break;
    case GxCases::QUADRATIC:
        get_quadratic_case(expected, g_diff);
        break;
    case GxCases::ABS_VALUE:
        get_abs_value_case(expected, g_diff);
        break;
    case GxCases::CIRCLE:
        get_circle_case(expected, g_diff);
        break;
    case GxCases::MAX:
        get_max_case(expected, g_diff);
        break;
    case GxCases::SMOOTH_MAX:
        get_smooth_max_case(expected, g_diff);
        break;
    }

    auto f = [&A, &b](const Eigen::VectorXd& x) -> double {
        return (A * x - b).squaredNorm() / 2;
    };

    auto grad_f = [&A, &b](const Eigen::VectorXd& x) -> Eigen::VectorXd {
        return (A * x - b);
    };

    auto g = [&g_diff](const Eigen::VectorXd x) -> Eigen::VectorXd {
        DVector gx = g_diff(x);
        Eigen::VectorXd g(gx.rows());
        for (int i = 0; i < gx.rows(); ++i) {
            g(i) = gx(i).getValue();
        }

        return g;
    };

    auto jac_g
        = [&g_diff, &g, &gx_case](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
        DVector gx = g_diff(x);
        Eigen::MatrixXd jac_gx(gx.rows(), NUM_VARS);
        for (int i = 0; i < gx.rows(); ++i) {
            jac_gx.row(i) = gx(i).getGradient();
        }
        Eigen::MatrixXd fd_jac_gx;
        ccd::finite_jacobian(x, g, fd_jac_gx);
        if (!ccd::compare_jacobian(jac_gx, fd_jac_gx)) {

            spdlog::warn("CASE {} CHECK_GRADIENT FAILED\n"
                         "\tgrad={}\n"
                         "\tfd  ={}",
                GxCasesNames[static_cast<int>(gx_case)], fmt_eigen(jac_gx),
                fmt_eigen(fd_jac_gx));
        }
        return jac_gx;
    };

    // setup callback
    std::vector<Eigen::VectorXd> it_x, it_alpha;
    std::vector<double> it_gamma;
    auto callback = [&it_x, &it_alpha, &it_gamma](const Eigen::VectorXd& x,
                        const Eigen::VectorXd& alpha, const double gamma) {
        it_x.push_back(x);
        it_alpha.push_back(alpha);
        it_gamma.push_back(gamma);
    };

    // make the call
    Eigen::VectorXd x(NUM_VARS), alpha(NUM_CONSTRAINTS);
    bool success = solve_ncp(A, b, g, jac_g, /*max_iter=*/300, callback,
        update_type, lcp_solver, x, alpha, check_convergence,
        /*check_convergence_unfeasible=*/false);

    // ------------------------------------------------------------------------
    // DEBUG INFO
    // ------------------------------------------------------------------------
    // --- save to file
    std::ofstream o(outfile);

    o << "it,x1,x2,alpha1,alpha2,gamma,f(xi),g(xi)1,g(xi)2" << std::endl;

    for (uint i = 0; i < it_x.size(); i++) {
        o << i << ",";
        o << it_x[i].format(CommaFmt) << ",";
        o << it_alpha[i].format(CommaFmt) << ",";
        o << it_gamma[i] << ",";
        o << f(it_x[i]) << ",";
        o << g(it_x[i]).format(CommaFmt) << ",";
        o << jac_g(it_x[i]).format(CommaFmt) << std::endl;
    }
    // expected
    o << "expected"
      << ",";
    o << expected.format(CommaFmt) << ",";
    o << "0.0,0.0"
      << ",";
    o << "0.0"
      << ",";
    o << f(expected) << ",";
    o << g(expected).format(CommaFmt) << ",";
    o << jac_g(expected).format(CommaFmt) << std::endl;

    // --- Log
    auto fx_is_optimal = std::abs(f(x) - f(expected)) < 1E-16;

    auto fx = f(x);
    auto fe = f(expected);
    std::stringstream ssx;
    ssx << std::setprecision(10) << x.format(CommaFmt);

    // Check KKT conditions
    // (1) grad f = grad g.T \alpha
    // (2) g(x) >= 0
    // (3) \alpha.T g(x) == 0
    Eigen::VectorXd kkt_grad = grad_f(x) - jac_g(x).transpose() * alpha;
    double kkt_lcp = alpha.transpose() * g(x);

    spdlog::info("CASE={} UPDATE_TYPE={} LCP_SOLVER={} \n"
                 "\tnum_it:{}\n"
                 "\tsucess:{}\n"
                 "\toptimal:{}\n"
                 "\tx*:{}\n"
                 "\te*:{}\n"
                 "\tf(x*):{:.10f}\n"
                 "\tf(e*):{:.10f}\n"
                 "\tg(x*):{}\n"
                 "\tg(e*):{}\n"
                 "\tkkt_grad:{}\n"
                 "\tkkt_lcp:{}",
        GxCasesNames[static_cast<int>(gx_case)],
        NcpUpdateNames[static_cast<int>(update_type)],
        LCPSolverNames[lcp_solver], it_x.size(), fmt_bool(success),
        fmt_bool(fx_is_optimal), fmt_eigen(x), fmt_eigen(expected), fx, fe,
        fmt_eigen(g(x)), fmt_eigen(g(expected)), fmt_eigen(kkt_grad), kkt_lcp);
}

int main(int argc, char* argv[])
{
    std::string out_dir = std::string(DATA_OUTPUT_DIR);
    if (argc >= 2) {
        out_dir = argv[1];
    }

    using namespace ccd::opt;
    LCPSolver lcp_solvers[] = { LCP_GAUSS_SEIDEL, LCP_MOSEK };
    NcpUpdate update_types[] = { NcpUpdate::G_GRADIENT, NcpUpdate::LINEARIZED };
    bool check_convergence = true;

    for (const auto& gx_case : GxCasesAll) {
        for (const auto& update_type : update_types) {
            for (const auto& lcp_solver : lcp_solvers) {
                std::string case_name = fmt::format("ncp_case_{0}_{1}_{2}.csv",
                    GxCasesNames[static_cast<int>(gx_case)],
                    NcpUpdateNames[static_cast<int>(update_type)],
                    LCPSolverNames[lcp_solver]);

                std::string filename = out_dir + "/" + case_name;
                eval_ncp(gx_case, update_type, lcp_solver, check_convergence,
                    filename);
            }
        }
    }
    return 0;
}
