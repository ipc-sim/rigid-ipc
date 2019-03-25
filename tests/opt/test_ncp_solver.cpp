#include <opt/ncp_solver.hpp>

#include <fstream>
#include <iomanip> // std::setw
#include <iostream>

#include <catch.hpp>
#include <fmt/format.h>
#include <spdlog/spdlog.h>

#include <autodiff/autodiff_types.hpp>
namespace ccd {
namespace unittests {

    std::string fmt_red(const std::string& s)
    {
        return "\033[1;31m" + s + "\033[0m";
    }

    std::string fmt_bool(const bool s)
    {
        if (!s) {
            return fmt_red("F");
        }
        return "T";
    }
    static const Eigen::IOFormat CommaFmt(Eigen::StreamPrecision,
        Eigen::DontAlignCols, ", ", ", ", "", "", "", "");

    std::string fmt_eigen(const Eigen::MatrixXd& x, const int precision = 10)
    {
        std::stringstream ssx;
        ssx << std::setprecision(precision) << x.format(CommaFmt);
        return ssx.str();
    }

    typedef std::function<double(const Eigen::VectorXd& x)> callback_f;
    typedef std::function<Eigen::VectorXd(const Eigen::VectorXd& x)> callback_g;

    double eval_f(const Eigen::SparseMatrix<double> A, const Eigen::VectorXd& b,
        const Eigen::VectorXd& x)
    {
        return (A * x - b).squaredNorm() / 2;
    }

    Eigen::VectorXd eval_grad_f(const Eigen::SparseMatrix<double> A,
        const Eigen::VectorXd& b, const Eigen::VectorXd& x)
    {
        return (A * x - b);
    }

    // differentiable helpers
    template <typename T> using VectorXT = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    template <typename T, int N> using VectorNT = Eigen::Matrix<T, N, 1>;
    template <int N> using DScalar = DScalar1<double, VectorNT<double, N>>;

    template <typename T> Eigen::VectorXd eval_g(const VectorXT<T>& gx)
    {
        Eigen::VectorXd g(gx.rows());
        for (uint i = 0; i < gx.rows(); ++i) {
            g(i) = gx(i).getValue();
        }

        return g;
    }

    template <typename T>
    Eigen::MatrixXd eval_jac_g(const VectorXT<T>& gx, const int num_vars)
    {
        Eigen::MatrixXd jac_gx(gx.rows(), num_vars);
        for (uint i = 0; i < gx.rows(); ++i) {
            jac_gx.row(i) = gx(i).getGradient();
        }

        return jac_gx;
    }

    // ---------------------------------------------------
    // Tests
    // ---------------------------------------------------
    TEST_CASE("NCP 0 gradient_only", "[opt][NCP][NCP-Interface][gradient_only]")
    {

        auto constraint = GENERATE(0, 1, 2, 3, 4);
        auto update_type = GENERATE(ccd::opt::GRADIENT_ONLY, ccd::opt::OTHER);
        auto lcp_solver
            = GENERATE(ccd::opt::LCP_GAUSS_SEIDEL, ccd::opt::LCP_MOSEK);

        using namespace ccd::opt;
        const int num_vars = 2;
        const int num_constraints = 2;

        Eigen::SparseMatrix<double> A(num_vars, num_vars);
        Eigen::VectorXd b(num_vars), expected(num_vars), x(num_vars),
            alpha(num_constraints);

        DiffScalarBase::setVariableCount(size_t(num_vars));
        typedef DScalar<num_vars> dscalar;

        // ------------------- Custom part
        // Min (x1 + 1)^2 + (x2 + 2.5)^2
        // s.t
        A.setIdentity();
        b << -1, -2.5;
        const char* tests[]
            = { "LINEAR", "QUADRATIC", "ABS_VALUE", "CIRCLE", "MAX" };

        auto g_diff
            = [constraint](const Eigen::VectorXd& x) -> VectorXT<dscalar> {
            VectorXT<dscalar> gx(num_constraints);
            dscalar x0(0, x[0]);
            dscalar x1(1, x[1]);

            switch (constraint) {
            default:
            case 0: // x0, x1 >= 0
                gx(0) = x0;
                gx(1) = x1;
                break;
            case 1: // x0^2 <= 0.04, x1^2 <= 0.09
                gx(0) = 0.04 - x0 * x0;
                gx(1) = 0.09 - x1 * x1;
                break;
            case 2: // |x0| <= 0.2, |x1| <= 0.3
                gx(0) = 0.2 - (x0 > 0 ? x0 : -x0);
                gx(1) = 0.3 - (x1 > 0 ? x1 : -x1);
                break;
            case 3: // (x0 - 1)^2 <= 1, (x1 - 2.5)^2 <= 1
                gx(0) = 1.0 - (x0 - 1.0) * (x0 - 1.0);
                gx(1) = 1.0 - (x1 - 2.5) * (x1 - 2.5);
                break;
            case 4: // max(x1^2,1) <=1, max(x2^2,1) <=1
                auto aux1 = x1 * x1;
                auto aux0 = x0 * x0;
                gx(0) = 1.0 - (aux0 > 1 ? aux0 : dscalar(1.0));
                gx(1) = 1.0 - (aux1 > 1 ? aux1 : dscalar(1.0));
                break;
            }
            return gx;
        };
        switch (constraint) {
        default:
        case 0:
            expected << 0.0, 0.0;
            break;
        case 1:
        case 2:
            expected << -0.2, -0.3;
            break;
        case 3:
            expected << 0.0, 1.5;
            break;
        case 4:
            expected << -1.0, -1.0;
            break;
        }

        // ---------------------- End Custom part

        auto f = [&A, &b](const Eigen::VectorXd& x) -> double {
            return eval_f(A, b, x);
        };

        auto grad_f = [&A, &b](const Eigen::VectorXd& x) -> Eigen::VectorXd {
            return eval_grad_f(A, b, x);
        };

        callback_g g = [g_diff](const Eigen::VectorXd x) -> Eigen::VectorXd {
            return eval_g(g_diff(x));
        };

        callback_jac_g jac_g
            = [g_diff](const Eigen::VectorXd& x) -> Eigen::MatrixXd {
            return eval_jac_g(g_diff(x), num_vars);
        };

        std::vector<Eigen::VectorXd> it_x;
        std::vector<Eigen::VectorXd> it_alpha;
        std::vector<double> it_gamma;
        callback_intermediate_ncp callback
            = [&it_x, &it_alpha, &it_gamma](const Eigen::VectorXd& x,
                  const Eigen::VectorXd& alpha, const double gamma) {
                  it_x.push_back(x);
                  it_alpha.push_back(alpha);
                  it_gamma.push_back(gamma);
              };

        bool success = solve_ncp(A, b, g, jac_g, /*max_iter=*/100, callback,
            update_type, lcp_solver, x, alpha);

        // ------------------------------------------------------------------------
        // DEBUG INFO
        // ------------------------------------------------------------------------
        // --- save to file
        const std::string base_dir = TEST_OUTPUT_DIR;
        std::string test_name = fmt::format("ncp_case_{0}_{1}_{2}.csv",
            tests[constraint], update_type, LCPSolverNames[lcp_solver]);
        std::ofstream o(base_dir + "/" + test_name);

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
        auto diff = (x - expected).squaredNorm();
        auto fx_is_optimal = f(x) == Approx(f(expected)).epsilon(1E-16);

        spdlog::info("CASE={} UPDATE_TYPE={} LCP_SOLVER={}\n"
                     "\tnum_it:{}\n"
                     "\tsucess:{}\n"
                     "\toptimal:{}",
            tests[constraint], update_type, LCPSolverNames[lcp_solver],
            it_x.size(), fmt_bool(success), fmt_bool(fx_is_optimal));

        auto fx = f(x);
        auto fe = f(expected);
        std::stringstream ssx;
        ssx << std::setprecision(10) << x.format(CommaFmt);

        spdlog::info("Actual vs Expected\n"
                     "\tx*:{}\n"
                     "\te*:{}\n"
                     "\tf(x*):{:.10f}\n"
                     "\tf(e*):{:.10f}",
            fmt_eigen(x), fmt_eigen(expected), fx, fe);

        // ------------------------------------------------------------------------
        // CHECK OPTIMALITY CONDITIONS
        // ------------------------------------------------------------------------
        REQUIRE(success);
        if (success && (diff > 1e-16 || !fx_is_optimal)) {
            // if results are different, the expected one should be smaller
            // otherwise expected value is not optimal.
            REQUIRE(fx > fe);
        }
        REQUIRE(fx_is_optimal);

        // Check KKT conditions
        // (1) grad f = grad g.T \alpha
        // (2) g(x) >= 0
        // (3) \alpha.T g(x) == 0
        Eigen::VectorXd kkt_grad = grad_f(x) - jac_g(x).transpose() * alpha;
        double kkt_cp = alpha.transpose() * g(x);

        CHECK((g(x).array() >= 0).all());
        CHECK(kkt_grad.squaredNorm() < 1e-16);
        CHECK(std::abs(kkt_cp) < 1E-16);
    }
} // namespace unittests
} // namespace ccd
