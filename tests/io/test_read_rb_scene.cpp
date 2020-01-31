#include <catch2/catch.hpp>

#include <iostream>

#include <io/read_rb_scene.hpp>
#include <physics/mass.hpp>

//-----------------
// Tests
//-----------------

TEST_CASE("Read 2D RB Scene", "[io][json][rigid-body]")
{
    std::vector<ccd::physics::RigidBody> rbs;
    using namespace nlohmann;
    auto j = R"({"rigid_bodies": [
             {
               "vertices": [
                 [0, 0],
                 [2, 0],
                 [2, 3],
                 [5, 3],
                 [5, 5],
                 [0, 5]
               ],
               "edges": [
                 [0, 1],
                 [1, 2],
                 [2, 3],
                 [3, 4],
                 [4, 5],
                 [5, 0]
               ],
               "linear_velocity": [0.0, 0.0],
               "angular_velocity": [0.0]
             }
           ]})"_json;

    ccd::io::read_rb_scene(j, rbs);

    CHECK(rbs.size() == 1);

    Eigen::MatrixXd positions(6, 2);
    positions << 0, 0, 2, 0, 2, 3, 5, 3, 5, 5, 0, 5;
    CHECK((positions - rbs[0].world_vertices()).squaredNorm() == 0.0);

    Eigen::MatrixXi edges(6, 2);
    edges << 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 0;
    CHECK((edges - rbs[0].edges).squaredNorm() == 0.0);

    Eigen::Vector3d velocity(3);
    velocity << 0, 0, 0;
    CHECK((velocity - rbs[0].velocity.dof()).squaredNorm() == 0.0);
}

// TODO: Test reading 3D RB scenes
