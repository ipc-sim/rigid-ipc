#include <catch2/catch.hpp>

#include <igl/PI.h>

#include <interval/interval.hpp>
#include <logger.hpp>
#include <physics/pose.hpp>

using namespace ipc;
using namespace ipc::rigid;

TEST_CASE("Simple interval arithmetic", "[interval]")
{
    ipc::rigid::Interval i(0, 1), j(4, 5);
    CHECK(i.lower() <= i.upper());
    CHECK(j.lower() <= j.upper());

    ipc::rigid::Interval result;

    result = i + j;
    CHECK(result.lower() <= result.upper());

    result = i - j;
    CHECK(result.lower() <= result.upper());

    result = i * j;
    CHECK(result.lower() <= result.upper());

    result = 1.0 / j;
    CHECK(result.lower() <= result.upper());

    result = i / j;
    CHECK(result.lower() <= result.upper());

    result = i + 10.0;
    CHECK(result.lower() == Approx(i.lower() + 10).margin(1e-12));
    CHECK(result.upper() == Approx(i.upper() + 10).margin(1e-12));
}

TEST_CASE("Cosine interval arithmetic", "[interval]")
{
    ipc::rigid::Interval r;

    double shift;
    SECTION("No shift") { shift = 0; }
    SECTION("2π shift") { shift = 2 * igl::PI; }
    SECTION("-2π shift") { shift = -2 * igl::PI; }
    SECTION("100π shift") { shift = 100 * igl::PI; }
    SECTION("-100π shift") { shift = -100 * igl::PI; }

    CAPTURE(shift);

    r = cos(ipc::rigid::Interval(-1, 7) + shift);
    CHECK(r.lower() == -1.0);
    CHECK(r.upper() == 1.0);

    r = cos(ipc::rigid::Interval(2, 4) + shift);
    CHECK(r.upper() < 0);

    r = cos(ipc::rigid::Interval(0, 1) + shift);
    CHECK(r.lower() > 0);

    r = cos(ipc::rigid::Interval(1, 2) + shift);
    CHECK(r.lower() < 0);
    CHECK(r.upper() > 0);
}

TEST_CASE("Sine interval arithmetic", "[interval]")
{
    ipc::rigid::Interval r;

    double shift = igl::PI / 2;
    SECTION("No shift") { shift += 0; }
    SECTION("2π shift") { shift += 2 * igl::PI; }
    SECTION("-2π shift") { shift += -2 * igl::PI; }
    SECTION("100π shift") { shift += 100 * igl::PI; }
    SECTION("-100π shift") { shift += -100 * igl::PI; }

    CAPTURE(shift);

    r = sin(ipc::rigid::Interval(-1, 7) + shift);
    CHECK(r.lower() == -1.0);
    CHECK(r.upper() == 1.0);

    r = sin(ipc::rigid::Interval(2, 4) + shift);
    CAPTURE(
        (ipc::rigid::Interval(2, 4) + shift).lower(),
        (ipc::rigid::Interval(2, 4) + shift).upper());
    CHECK(r.upper() < 0);

    r = sin(ipc::rigid::Interval(0, 1) + shift);
    CHECK(r.lower() > 0);

    r = sin(ipc::rigid::Interval(1, 2) + shift);
    CHECK(r.lower() < 0);
    CHECK(r.upper() > 0);
}

TEST_CASE("Interval rotation rounding", "[interval][matrix]")
{
    Pose<double> pose(
        /*x=*/0.30969396267858817, /*y=*/0.85675409103416755,
        /*theta=*/0.79358805865013693);
    Eigen::Matrix<double, 4, 2> V;
    V.row(0) << 0.5, 0.0;
    V.row(1) << 0.0, 0.5;
    V.row(2) << -0.5, 0.0;
    V.row(3) << 0.0, -0.5;

    Pose<Interval> posei = pose.cast<Interval>();
    MatrixMax3<Interval> R = posei.construct_rotation_matrix();
    MatrixMax3<Interval> expected_R;
    expected_R.resize(2, 2);
    expected_R.row(0) << cos(posei.rotation(0)), -sin(posei.rotation(0));
    expected_R.row(1) << sin(posei.rotation(0)), cos(posei.rotation(0));

    CHECK((R - expected_R).squaredNorm().lower() == Approx(0.0).margin(1e-12));
    CHECK((R - expected_R).squaredNorm().upper() == Approx(0.0).margin(1e-12));

    Interval interval = R(0, 0);
    CHECK(!boost::numeric::empty(interval));
    Interval tmp = 0.5 * interval;

    MatrixMax3<Interval> RT = R.transpose();
    // Eigen::Matrix<Interval, 4, 2> RV = Vi * RT;
    Eigen::Matrix<Interval, 4, 2> RV;
    for (int i = 0; i < V.rows(); i++) {
        for (int j = 0; j < V.cols(); j++) {
            Interval tmp1 = V(i, j) * R(0, j);
            Interval tmp2 = V(i, j) * R(1, j);
            RV(i, 0) = tmp1 + tmp2;
        }
    }
    RV.rowwise() += posei.position.transpose();
    for (int i = 0; i < RV.size(); i++) {
        CHECK(std::isfinite(RV(i).lower()));
        CHECK(std::isfinite(RV(i).upper()));
    }
    // std::cout << RV << std::endl;
}

#ifdef USE_FILIB_INTERVALS
TEST_CASE("CheckingPolicy", "[interval][checking_policy]")
{
    typedef boost::numeric::interval_lib::checking_no_empty<double>
        m_CheckingPolicy;

    // Use filib rounding arithmetic
    typedef boost::numeric::interval<
        double,
        boost::numeric::interval_lib::policies<
            boost::numeric::interval_lib::save_state<FILibRounding>,
            m_CheckingPolicy>>
        m_Interval;

    m_Interval i = m_Interval(-1, 1);
    // fmt::print("{}\n", logger::fmt_interval(i / 0.0));
}
#endif
