#ifndef COLLISION_GENERATOR_HPP
#define COLLISION_GENERATOR_HPP

#include <memory>

#include <Eigen/Core>
#include <catch2/catch.hpp>

namespace ccd {
namespace unittests {

    struct TestImpact {
        Eigen::Vector2d Vi, Vj, Vk, Vl;
        Eigen::Vector2d Ui, Uj, Uk, Ul;
        double toi, alpha;
    };
    std::ostream& operator << (std::ostream& o, const TestImpact& p);

    TestImpact generate_random_impacts(const bool rigid);

    class TestImpactGenerator
        : public Catch::Generators::IGenerator<TestImpact> {
        bool rigid;

    public:
        TestImpactGenerator(bool rigid);

        auto get(size_t index) const -> TestImpact override;
    };

    Catch::Generators::Generator<TestImpact> random_impacts(
        size_t value, bool rigid);

}

}
#endif // COLLISION_GENERATOR_H
