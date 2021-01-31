#include <catch2/catch.hpp>

#include <solvers/line_search.cpp>

#include <logger.hpp>

TEST_CASE("Test unconstrained line search", "[opt][line_search]")
{
    auto f = [](const Eigen::VectorXd& x) -> double {
        return x.squaredNorm() / 2;
    };
    bool use_wolfe = GENERATE(false, true);
    auto grad_f = [](const Eigen::VectorXd& x) -> Eigen::VectorXd { return x; };
    Eigen::VectorXd x, dir;
    double expected_step_length = -1, actual_step_length = 1;
    bool expected_success = true, actual_success;
    SECTION("Half step")
    {
        x = (Eigen::VectorXd(1) << -1).finished();
        dir = (Eigen::VectorXd(1) << 2).finished();
        expected_step_length = 0.5;
    }
    SECTION("Full step")
    {
        x = (Eigen::VectorXd(1) << -1).finished();
        dir = (Eigen::VectorXd(1) << 1).finished();
        expected_step_length = 1.0;
    }
    SECTION("No step")
    {
        x = (Eigen::VectorXd(1) << -1).finished();
        dir = (Eigen::VectorXd(1) << -2).finished();
        expected_success = false;
    }
    if (use_wolfe) {
        actual_success
            = ccd::opt::line_search(x, dir, f, grad_f(x), actual_step_length);
    } else {
        actual_success = ccd::opt::line_search(x, dir, f, actual_step_length);
    }
    CHECK(actual_success == expected_success);
    if (expected_success) {
        CHECK(actual_step_length == Approx(expected_step_length));
    }
}

TEST_CASE("Test constrained line search", "[opt][line_search]")
{
    auto f = [](const Eigen::VectorXd& x) -> double {
        return x.squaredNorm() / 2;
    };
    auto grad_f = [](const Eigen::VectorXd& x) -> Eigen::VectorXd { return x; };
    auto constraint
        = [](const Eigen::VectorXd& x) -> bool { return x(0) >= 0.5; };

    Eigen::VectorXd x, dir;
    double expected_step_length = -1, actual_step_length = 1;
    bool expected_success = true, actual_success;
    SECTION("Half step")
    {
        x = (Eigen::VectorXd(1) << 1).finished();
        dir = (Eigen::VectorXd(1) << -1).finished();
        expected_step_length = 0.5;
    }
    SECTION("Full step")
    {
        x = (Eigen::VectorXd(1) << 1).finished();
        dir = (Eigen::VectorXd(1) << -0.5).finished();
        expected_step_length = 1.0;
    }
    SECTION("No step")
    {
        x = (Eigen::VectorXd(1) << 0.5).finished();
        dir = (Eigen::VectorXd(1) << -0.5).finished();
        expected_success = false;
    }
    actual_success = ccd::opt::constrained_line_search(
        x, dir, f, grad_f(x), constraint, actual_step_length);
    CAPTURE(actual_step_length);
    CHECK(actual_success == expected_success);
    if (expected_success) {
        CHECK(actual_step_length == Approx(expected_step_length));
    }
}
