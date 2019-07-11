#include <array>
#include <iomanip>
#include <iostream>

#include <catch2/catch.hpp>

#include <igl/PI.h>

#include <physics/rigid_body_assembler.hpp>

// ---------------------------------------------------
// Tests
// ---------------------------------------------------

TEST_CASE("Rigid Body System Transform", "[RB][RB-System][RB-System-transform]")
{

    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(4, 2);
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    Eigen::Vector3d rb_displ_1, rb_displ_2;

    Eigen::MatrixX2d expected(4, 2);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        rb_displ_1 << 0.5, 0.5, 0.0;
        rb_displ_2 << 1.0, 1.0, 0.0;

        // expected displacements
        expected.resize(8, 2);
        expected.block(0, 0, 4, 2)
            = rb_displ_1.segment(0, 2).transpose().replicate(4, 1);
        expected.block(4, 0, 4, 2)
            = rb_displ_2.segment(0, 2).transpose().replicate(4, 1);
    }

    SECTION("90 Deg Rotation Case")
    {
        rb_displ_1 << 0.0, 0.0, 0.5 * M_PI;
        rb_displ_2 << 0.0, 0.0, M_PI;
        expected.resize(8, 2);
        expected.block(0, 0, 4, 2) << 1.0, 0.0, 0.0, 1.0, -1.0, 0.0, 0.0, -1.0;
        expected.block(4, 0, 4, 2) << 1.0, 1.0, -1.0, 1.0, -1.0, -1.0, 1.0,
            -1.0;
    }

    using namespace ccd::physics;

    std::vector<RigidBody> rbs;
    RigidBodyAssembler assembler;
    rbs.push_back(RigidBody::from_velocity(vertices, edges, velocity));
    rbs.push_back(RigidBody::from_velocity(vertices, edges, velocity));
    assembler.init(rbs);

    Eigen::VectorXd pos(6);
    pos.segment(0, 3) = rb_displ_1 + rbs[0].position;
    pos.segment(3, 3) = rb_displ_2 + rbs[1].position;

    /// compute displacements between current and given positions
    /// TODO: update test to not need displacements
    Eigen::MatrixXd actual = assembler.world_vertices(pos) - assembler.world_vertices();
    CHECK((expected - actual).squaredNorm() < 1E-6);

    //    Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision,
    //    Eigen::DontAlignCols,
    //        ", ", ", ", "", "", " << ", ";");
    //    std::cout << expected.format(CommaInitFmt) << std::endl;
    //    std::cout << actual.format(CommaInitFmt) << std::endl;
}

TEST_CASE("Rigid Body System Gradient", "[RB][RB-System][RB-System-gradient]")
{

    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(4, 2);
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    Eigen::Vector3d rb_displ_1, rb_displ_2;

    Eigen::MatrixXd expected(16, 6);

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        rb_displ_1 << 0.5, 0.5, 0.0;
        rb_displ_2 << 1.0, 1.0, 0.0;

        // expected displacements
        for (uint i = 0; i < 4; ++i) {
            // x - derivative RB 1
            expected.row(i) << 1, 0, -vertices(i, 1), 0, 0, 0;
            // x - derivative RB 2
            expected.row(4 + i) << 0, 0, 0, 1, 0, -vertices(i, 1);
            // y - derivative RB 1
            expected.row(8 + i) << 0, 1, vertices(i, 0), 0, 0, 0;
            // y - derivative RB 1
            expected.row(12 + i) << 0, 0, 0, 0, 1, vertices(i, 0);
        }
    }

    SECTION("90 Deg Rotation Case")
    {
        rb_displ_1 << 0.0, 0.0, 0.5 * M_PI;
        rb_displ_2 << 0.0, 0.0, M_PI;

        // expected displacements
        for (uint i = 0; i < 4; ++i) {
            // x - derivative RB 1
            double dtheta = -vertices(i, 0) * sin(rb_displ_1(2))
                - vertices(i, 1) * cos(rb_displ_1(2));
            expected.row(i) << 1, 0, dtheta, 0, 0, 0;
            // x - derivatibe RB 2
            dtheta = -vertices(i, 0) * sin(rb_displ_2(2))
                - vertices(i, 1) * cos(rb_displ_2(2));
            expected.row(4 + i) << 0, 0, 0, 1, 0, dtheta;
            // y - derivatibe RB 1
            dtheta = vertices(i, 0) * cos(rb_displ_1(2))
                - vertices(i, 1) * sin(rb_displ_1(2));
            expected.row(8 + i) << 0, 1, dtheta, 0, 0, 0;
            // y - derivative RB 1
            dtheta = vertices(i, 0) * cos(rb_displ_2(2))
                - vertices(i, 1) * sin(rb_displ_2(2));
            expected.row(12 + i) << 0, 0, 0, 0, 1, dtheta;
        }
    }

    using namespace ccd::physics;


    std::vector<RigidBody> rbs;
    RigidBodyAssembler assembler;
    rbs.push_back(RigidBody::from_velocity(vertices, edges, velocity));
    rbs.push_back(RigidBody::from_velocity(vertices, edges, velocity));
    assembler.init(rbs);

    Eigen::VectorXd pos(6);
    pos.segment(0, 3) = rb_displ_1 + rbs[0].position;
    pos.segment(3, 3) = rb_displ_2 + rbs[1].position;

    Eigen::SparseMatrix<double> actual;
    assembler.world_vertices_gradient(pos, actual);
    CHECK((expected - actual.toDense()).squaredNorm() < 1E-6);

    //    Eigen::IOFormat HeavyFmt(4, 0, ", ", ";\n", "[", "]", "[", "]");
    //    std::cout << expected.format(HeavyFmt) << std::endl;
    //    std::cout << actual.toDense().format(HeavyFmt) << std::endl;
}

TEST_CASE("Rigid Body System Hessian", "[RB][RB-System][RB-System-hessian]")
{

    Eigen::MatrixX2d vertices(4, 2);
    Eigen::MatrixX2i edges(4, 2);
    Eigen::Vector3d velocity = Eigen::Vector3d::Zero();
    Eigen::Vector3d rb_displ_1, rb_displ_2;

    std::array<Eigen::Matrix<double, 6, 6>, 16> expected;

    vertices << -0.5, -0.5, 0.5, -0.5, 0.5, 0.5, -0.5, 0.5;
    edges << 0, 1, 1, 2, 2, 3, 3, 0;

    SECTION("Translation Case")
    {
        rb_displ_1 << 0.5, 0.5, 0.0;
        rb_displ_2 << 1.0, 1.0, 0.0;
    }
    SECTION("90 Deg Rotation Case")
    {
        rb_displ_1 << 0.0, 0.0, 0.5 * M_PI;
        rb_displ_2 << 0.0, 0.0, M_PI;
    }

    // clang-format off
        for (size_t i = 0; i < 4; ++i) {
            // x body 1
            expected[i].setZero();
            expected[i](2,2) = - vertices(i,0) * cos(rb_displ_1(2)) + vertices(i,1) * sin(rb_displ_1(2));

            // y body 1
            expected[8 + i].setZero();
            expected[8 + i](2,2) = - vertices(i,0) * sin(rb_displ_1(2)) - vertices(i,1) * cos(rb_displ_1(2));

            // x body 2
            expected[4 + i].setZero();
            expected[4 + i](5,5) = - vertices(i,0) * cos(rb_displ_2(2)) + vertices(i,1) * sin(rb_displ_2(2));

            // y body 2
            expected[12 + i].setZero();
            expected[12 + i](5,5) = - vertices(i,0) * sin(rb_displ_2(2)) - vertices(i,1) * cos(rb_displ_2(2));

        }

    // clang-format on

    using namespace ccd::physics;

    std::vector<RigidBody> rbs;
    RigidBodyAssembler assembler;
    rbs.push_back(RigidBody::from_velocity(vertices, edges, velocity));
    rbs.push_back(RigidBody::from_velocity(vertices, edges, velocity));
    assembler.init(rbs);

    Eigen::VectorXd pos(6);
    pos.segment(0, 3) = rb_displ_1 + rbs[0].position;
    pos.segment(3, 3) = rb_displ_2 + rbs[1].position;

    std::vector<Eigen::SparseMatrix<double>> actual;
    assembler.world_vertices_hessian(pos, actual);
    REQUIRE(actual.size() == expected.size());
    for (size_t i = 0; i < actual.size(); ++i) {
        REQUIRE(expected[i].rows() == actual[i].rows());
        REQUIRE(expected[i].cols() == actual[i].cols());
//        std::cout << i << std::endl;
//        Eigen::IOFormat HeavyFmt(4, 0, ", ", ";\n", "[", "]", "[", "]");
//        std::cout << expected[i].format(HeavyFmt) << std::endl;
//        std::cout << actual[i].toDense().format(HeavyFmt) << std::endl;
        REQUIRE((expected[i] - actual[i].toDense()).squaredNorm() < 1E-6);
    }
}
