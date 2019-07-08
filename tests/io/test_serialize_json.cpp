#include <catch2/catch.hpp>

#include <io/serialize_json.hpp>
#include <iostream>

//-----------------
// Tests
//-----------------

TEST_CASE("Serialize to Json", "[io][json]")
{

    SECTION("Serialize double vector")
    {
        Eigen::VectorXd vec(3);
        vec << 1.0, 2.0, 3.0;
        auto jvec = ccd::io::to_json<double>(vec);

        CHECK(jvec.dump() == "[1.0,2.0,3.0]");
    }

    SECTION("Serialize int vector")
    {
        Eigen::VectorXi vec(3);
        vec << 1, 2, 3;
        auto jvec = ccd::io::to_json<int>(vec);

        CHECK(jvec.dump() == "[1,2,3]");
    }

    SECTION("Serialize double Row Vector")
    {
        Eigen::MatrixXd vec(1, 3);
        vec << 1.0, 2.0, 3.0;
        auto jvec = ccd::io::to_json<double>(vec);

        CHECK(jvec.dump() == "[[1.0,2.0,3.0]]");
    }

    SECTION("Serialize double matrix")
    {
        Eigen::MatrixXd vec(2, 3);
        vec << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
        auto jvec = ccd::io::to_json<double>(vec);

        CHECK(jvec.dump() == "[[1.0,2.0,3.0],[4.0,5.0,6.0]]");
    }

    SECTION("Serialize int Row Vector")
    {
        Eigen::MatrixXi vec(1, 3);
        vec << 1, 2, 3;
        auto jvec = ccd::io::to_json<int>(vec);

        CHECK(jvec.dump() == "[[1,2,3]]");
    }
}

TEST_CASE("Deserialize from Json", "[io][json]")
{

    SECTION("Deserialize double vector")
    {
        Eigen::VectorXd actual(3), expected(3);
        expected << 1.0, 2.0, 3.0;
        nlohmann::json json = "[1.0, 2.0, 3.0]"_json;
        ccd::io::from_json<double>(json, actual);

        CHECK((actual - expected).squaredNorm() == 0.0);
    }

    SECTION("Deserialize int vector")
    {
        Eigen::VectorXi actual(3), expected(3);
        expected << 1, 2, 3;
        nlohmann::json json = "[1, 2, 3]"_json;
        ccd::io::from_json<int>(json, actual);

        CHECK((actual - expected).squaredNorm() == 0.0);
    }

    SECTION("Deserialize double Row Vector")
    {
        Eigen::MatrixXd actual(1, 3), expected(1, 3);
        expected << 1.0, 2.0, 3.0;
        nlohmann::json json = "[[1.0,2.0,3.0]]"_json;
        ccd::io::from_json<double>(json, actual);

        CHECK((actual - expected).squaredNorm() == 0.0);
    }

    SECTION("Deserialize double matrix")
    {
        Eigen::MatrixXd actual(2, 3), expected(2, 3);
        expected << 1.0, 2.0, 3.0, 4.0, 5.0, 6.0;
        nlohmann::json json = "[[1.0,2.0,3.0],[4.0,5.0,6.0]]"_json;
        ccd::io::from_json<double>(json, actual);

        CHECK((actual - expected).squaredNorm() == 0.0);
    }

    SECTION("Deserialize int Row Vector")
    {
        Eigen::MatrixXi actual(1, 3), expected(1, 3);
        expected << 1, 2, 3;
        nlohmann::json json = "[[1,2,3]]"_json;
        ccd::io::from_json<int>(json, actual);

        CHECK((actual - expected).squaredNorm() == 0.0);
    }
}
