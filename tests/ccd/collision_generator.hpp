#ifndef COLLISION_GENERATOR_HPP
#define COLLISION_GENERATOR_HPP

#include <memory>

#include <Eigen/Core>
#include <catch2/catch.hpp>

namespace ipc::rigid {
namespace unittests {

    struct TestImpact {
        Eigen::Vector2d Vi, Vj, Vk, Vl;
        Eigen::Vector2d Ui, Uj, Uk, Ul;
        double toi, alpha;
    };
    std::ostream& operator<<(std::ostream& o, const TestImpact& p);

    TestImpact generate_random_impact(const bool rigid);

    class TestImpactGenerator
        : public Catch::Generators::IGenerator<TestImpact> {

        TestImpact current;
        bool rigid;
        size_t current_i, max_i;

    public:
        TestImpactGenerator(size_t value, bool rigid);

        bool next() override;

        TestImpact const& get() const override;
    };

    Catch::Generators::GeneratorWrapper<TestImpact>
    random_impacts(size_t value, bool rigid);

} // namespace unittests

} // namespace ipc::rigid
#endif // COLLISION_GENERATOR_H
